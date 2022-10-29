import logging
import pdb

import ee

from sla.settings import ALOS, ING, SRTM

logger = logging.getLogger(__name__)


def extract_sla_patch(image):
    """
    Extracts Snow Line Altitute (SLA) patch from an image.
    """

    def calculate_area(img):
        """
        Helper method to calculate the area of a band without masked pixels.
        """

        # Take a band (elevation in this case) to calculate the area
        pixel_area = img.select("AVE_DSM").multiply(ee.Image.pixelArea())
        return pixel_area.reduceRegion(
            **{
                "reducer": ee.Reducer.sum(),  # type: ignore
                "geometry": geometry,
                "scale": 20,
                "maxPixels": 1e8,
                "bestEffort": True,
            }
        ).getNumber("AVE_DSM")

    ing_id = image.get("ING_ID")
    otsu = image.get("otsu")
    dem_info = image.get("deminfo")

    geometry = ING.filterMetadata("ING_ID", "equals", ing_id)

    # DEM Selection
    dem_1 = SRTM.select("elevation").rename("AVE_DSM").clip(geometry)
    dem_2 = ALOS.select("AVE_DSM").clip(geometry)

    dem_selector = ee.Algorithms.IsEqual(ee.String(dem_info), ee.String("ALOS"))
    dem_glacier = ee.Image(ee.Algorithms.If(dem_selector, dem_2, dem_1))

    # Combine classes from classified image

    classified = image
    dcmask = classified.select("classification").neq(-1)  # Unclassified
    mask1 = classified.select("classification").neq(7)  # All but shadow on water
    mask2 = classified.select("classification").neq(4)  # All but clouds
    mask3 = classified.select("classification").neq(3)  # All but debris cover
    mask4 = classified.select("classification").neq(2)  # All but water
    mask5 = classified.select("classification").neq(8)

    classified = (
        classified.mask(dcmask)
        .clip(geometry)
        .updateMask(mask1)
        .clip(geometry)
        .updateMask(mask2)
        .clip(geometry)
        .updateMask(mask3)
        .clip(geometry)
        .updateMask(mask4)
        .clip(geometry)
        .updateMask(mask5)
        .clip(geometry)
    )
    # TODO: EXPORT
    # Elevation analysis

    elevation = dem_glacier.select("AVE_DSM").clip(geometry)
    # TODO: EXPORT
    # Combine snow and shadow on snow (sos)
    single_class_snow = classified.select("classification").eq(1)
    single_class_sos = classified.select("classification").eq(6)
    snow_mask = single_class_snow.max(single_class_sos)

    # Combine ice and shadow on ice (soi)
    single_class_ice = classified.select("classification").eq(0)
    single_class_soi = classified.select("classification").eq(5)
    ice_mask = single_class_ice.max(single_class_soi)

    # Mask ice
    elev_class_ice = elevation.mask(ice_mask.clip(geometry))

    # Mask snow
    elev_class_snow = elevation.mask(snow_mask.clip(geometry))

    # Create raster image to vectorize
    snow_vec = (
        classified.mask(snow_mask)
        .select("classification")
        .multiply(0)
        .unmask(-10)
        .clip(geometry)
    )

    ice_vec = (
        classified.mask(ice_mask)
        .select("classification")
        .add(10)
        .unmask(-10)
        .clip(geometry)
    )

    vector_image = snow_vec.max(ice_vec)

    # Calculate and store Areas for Ratio calculations

    ice_area = calculate_area(elev_class_ice)
    # TODO: EXPORT
    snow_area = calculate_area(elev_class_snow)
    # TODO: EXPORT
    ice_snow_area = ice_area.add(snow_area)
    # TODO: EXPORT

    snow_part = snow_area.divide(ice_snow_area)
    void_part = ee.Number(1).subtract(
        (
            ee.Number(ice_area)
            .add(ee.Number(snow_area))
            .divide(calculate_area(elevation))
        )
    )

    # Set Snow condition to true when 95% of the snow/ice area is covered with snow
    # With this condition the SnowLine is set to the minimal altitude of the glacier-
    # outline area
    snow_cond = snow_part.gt(0.95)

    # Calculate lowest / highest possible SLA's

    lowest_sla = ee.Number(
        elevation.reduceRegion(
            **{
                "reducer": ee.Reducer.min(),  # type: ignore
                "geometry": geometry,
                "scale": 30,
                "bestEffort": True,
            }
        ).get("AVE_DSM")
    )

    highest_sla = ee.Number(
        elevation.reduceRegion(
            **{
                "reducer": ee.Reducer.max(),  # type: ignore
                "geometry": geometry,
                "scale": 30,
                "bestEffort": True,
            }
        ).get("AVE_DSM")
    )

    # Vectorize the classified map with snow and ice patches
    classes = vector_image.reduceToVectors(
        **{
            "reducer": ee.Reducer.countEvery(),  # type: ignore
            "geometry": geometry,
            "scale": 20,
            "eightConnected": False,
            "bestEffort": True,
            "maxPixels": 1e9,
        }
    )

    # Calculate the biggest snow ice patch
    snow_ice_vector_map = ee.FeatureCollection(classes)
    snow_filter = ee.Filter.eq("label", 0)
    snow_max_collection = snow_ice_vector_map.filter(snow_filter).sort("count", False)
    snow_area_vec = ee.Number(snow_max_collection.aggregate_sum("count"))
    snow_max = ee.Feature(snow_max_collection.first())

    ice_filter = ee.Filter.greaterThanOrEquals("label", 9)
    ice_max_collection = snow_ice_vector_map.filter(ice_filter).sort("count", False)
    ice_area_vec = ee.Number(
        ice_max_collection.aggregate_sum("count")
    )  # ice_are_vec is UNUSED as per lint
    # TODO: EXPORT
    ice_max = ee.Feature(ice_max_collection.first())
    # TODO: EXPORT

    # Check this implementation Essentially nulls are being replace with a conditional
    is_ice_null = ee.Feature(None).set("count", 0)
    ice_max = ee.Feature(ee.Algorithms.If(ice_max, ice_max, is_ice_null))

    snow_ice_fc = ee.Algorithms.Collection([snow_max, ice_max])
    snow_patch_area = ee.Algorithms.If(
        snow_max_collection.size(), ee.Number(snow_max.get("count")), 0
    )

    ice_patch_area = ee.Number(ice_max.get("count"))
    # TODO: EXPORT
    total_area = ee.Number(snow_ice_vector_map.aggregate_sum("count"))

    snow_area_ratio = snow_area_vec.divide(total_area)

    snow_patch_ratio = ee.Number(snow_patch_area).divide(total_area).multiply(100)
    ice_patch_ratio = ee.Number(ice_patch_area).divide(total_area).multiply(100)
    # TODO: EXPORT

    relevant_area = snow_patch_ratio.add(ice_patch_ratio)

    # Extract the zone where the two patches touch each other

    snow_ice_image = snow_ice_fc.reduceToImage(["label"], ee.Reducer.first())  # type: ignore
    # TODO: EXPORT

    bigger_ice = snow_ice_image.mask(snow_ice_image.select("first").eq(0)).focal_max(2)
    touching_zone = bigger_ice.subtract(
        snow_ice_image.mask(snow_ice_image.select("first").gte(9))
    )
    # TODO: EXPORT

    elev_touch = elevation.addBands(touching_zone, ["first"])
    elevation_snow_line = elev_touch.mask(elev_touch.select("first").lt(-1))

    # Calculate mean altitude of the zone where the patches touch
    mean_altitudes = elevation_snow_line.reduceRegion(
        **{
            "reducer": ee.Reducer.median(),  # type: ignore
            "maxPixels": 1e8,
            "geometry": geometry,
            "bestEffort": True,
        }
    )

    # Exception Handling

    lowest_snow_elev = ee.Number(
        elevation.clip(snow_max).reduceRegion(
            **{
                "reducer": ee.Reducer.min(),  # type: ignore
                "geometry": geometry,
                "scale": 30,
                "bestEffort": True,
            }
        )
    )

    no_touch = ee.Number(
        ee.Algorithms.If(snow_max_collection.size(), lowest_snow_elev, highest_sla)
    )

    snow_line_1 = mean_altitudes.get("AVE_DSM")
    snow_line_2 = ee.Algorithms.If(snow_line_1, snow_line_1, no_touch)

    snow_line = ee.Number(ee.Algorithms.If(snow_cond, lowest_sla, snow_line_2))

    snow_line_std_dev = ee.Number(
        elevation_snow_line.reduceRegion(
            **{
                "reducer": ee.Reducer.stdDev(),  # type: ignore
                "maxPixels": 1e8,
                "geometry": geometry,
                "bestEffort": True,
            }
        ).get("AVE_DSM")
    )

    # Check if 95% of the glacier is convered with snow by settings SLA to the lowest
    # elevation of the glacier
    coverage_95 = ee.Algorithms.If(snow_cond, 0, snow_line_std_dev)
    std_dev_sla = ee.Algorithms.If(mean_altitudes.size(), coverage_95, 0)

    # Add date in MS-Excel readable format
    image_date = ee.Number(classified.get("system_time_start"))

    shorten = image_date.divide(1000).floor().divide(86400).add(25569)

    # Create feature to return and export information in a table (eg CSV)
    feature = ee.Feature(None)
    feature = (
        feature.set("Snow Cover Ratio", snow_area_ratio)
        .set("SLA MP-Approach", snow_line)
        .set("ice_area", ice_area)
        .set("total_area", total_area)
        .set("Ratio Area v/o Snow or Ice", void_part)
        .set("system:time_start", image_date)
        .set("date_MSxlsx", shorten)
        .set("otsu", otsu)
        .set("Considered Area for MP-SLA Extraction", relevant_area)
        .set("stdDev SLA, 0 means no touch between Main patches [m]", std_dev_sla)
        .set("ID_Glacier", ing_id)
    )

    return feature
