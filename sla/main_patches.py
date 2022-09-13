import ee

RGI = ee.FeatureCollection("users/josiaszeller/RGVI_v6")
ALOS = ee.Image("JAXA/ALOS/AW3D30/V2_2")
SRTM = ee.Image("USGS/SRTMGL1_003")


def extract_sla_patch(image):
    """
    Extracts Snow Line Altitute (SLA) patch from an image.
    """

    def calculate_area(img):
        """
        Helper method to calculate the area of a band without masked pixels.
        """

        # Take a band (elevation in this case) to calculate the area
        prep = img.select("AVE_DSM")
        pixel_area = prep.multiply(ee.Image.pixelArea())
        pixel_area_reduced = pixel_area.reduceRegion(
            **{
                "reducer": ee.Reducer.sum(),  # type: ignore
                "geometry": geometry,
                "scale": 20,
                "maxPixels": 1e8,
                "bestEffort": True,
            }
        )

        return pixel_area_reduced.getNumber("AVE_DSM")

    glims_id = image.get("GLIMS_ID")
    otsu = image.get("otsu")
    dem_info = image.get("deminfo")

    # TODO CHECK IF THIS (SHOULD BE ING)
    geometry = RGI.filterMetadata("GLIMSId", "equals", glims_id)

    # DEM Selection
    dem_1 = SRTM.select("elevation").rename("AVE_DSM").clip(geometry)
    dem_2 = ALOS.select("AVE_DSM").clip(geometry)

    dem_selector = ee.Algorithms.IsEqual(ee.String(dem_info), ee.String("ALOS"))
    dem_glacier = ee.Image(ee.Algorithms.If(dem_selector, dem_2, dem_1))
    classified = image

    # Combine classes from classified image

    dcmask = classified.select("classification").neq(-1)  # Unclassified
    mask1 = classified.select("classification").neq(7)  # All but shadow on water
    mask2 = classified.select("classification").neq(4)  # All but clouds
    mask3 = classified.select("classification").neq(3)  # All but debris cover
    mask4 = classified.select("classification").neq(2)  # All but water
    mask5 = classified.select("classification").neq(8)

    classmask1 = classified.mask(dcmask).clip(geometry)
    classmask2 = classmask1.updateMask(mask1).clip(geometry)
    classmask3 = classmask2.updateMask(mask2).clip(geometry)
    classmask4 = classmask3.updateMask(mask3).clip(geometry)
    classmask5 = classmask4.updateMask(mask4).clip(geometry)
    classmask6 = classmask5.updateMask(mask5).clip(geometry)

    # Elevation analysis

    elevation = dem_glacier.select("AVE_DSM").clip(geometry)

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
    vector_image = (
        classified.mask(snow_mask)
        .select("classification")
        .multiply(0)
        .unmask(-10)
        .clip(geometry)
        .max(
            classified.mask(ice_mask)
            .select("classification")
            .add(10)
            .unmask(-10)
            .clip(geometry)
        )
    )

    # Calculate and store Areas for Ratio calculations

    ice_area = calculate_area(elev_class_ice)
    snow_area = calculate_area(elev_class_snow)
    ice_snow_area = ice_area.add(snow_area)

    snow_part = snow_area.divide(ice_snow_area)
    void_part = ee.Number(1).subtract(
        (
            ee.Number(ice_area)
            .add(ee.Number(snow_area))
            .divide(calculate_area(elevation))
        )
    )

    snow_ice_ratio = snow_area.divide(ice_area)
    ice_snow_ratio = ice_area.divide(snow_area)

    snow_cond = snow_part.gt(0.95)

    # Calculate lowest / highest possible SLA's

    lowest_sla = elevation.reduceRegion(
        **{
            "reducer": ee.Reducer.min(),  # type: ignore
            "geometry": geometry,
            "scale": 30,
            "bestEffort": True,
        }
    ).get("AVE_DSM")

    highest_sla = elevation.reduceRegion(
        **{
            "reducer": ee.Reducer.max(),  # type: ignore
            "geometry": geometry,
            "scale": 30,
            "bestEffort": True,
        }
    ).get("AVE_DSM")

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
    ice_area_vec = ee.Number(ice_max_collection.aggregate_sum("count"))
    ice_max = ee.Feature(ice_max_collection.first())

    # Check this implementation Essentially nulls are being replace with a conditional
    is_ice_null = ee.Feature(None).set("count", 0)
    ice_max = ee.Feature(ee.Algorithms.If(ice_max, ice_max, is_ice_null))

    snow_ice_fc = ee.Algorithms.Collection([snow_max, ice_max])
    snow_patch_area = ee.Algorithms.If(
        snow_max_collection.size(), ee.Number(snow_max.get("count")), 0
    )

    ice_patch_area = ee.Number(ice_max.get("count"))
    total_area = ee.Number(snow_ice_vector_map.aggregate_sum("count"))

    snow_area_ratio = snow_area_vec.divide(total_area)
    ice_area_ratio = ice_area_vec.divide(total_area)

    snow_patch_ratio = ee.Number(snow_patch_area).divide(total_area).multiply(100)
    ice_patch_ratio = ee.Number(ice_patch_area).divide(total_area).multiply(100)

    relevant_area = snow_patch_ratio.add(ice_patch_ratio)

    # Extract the zone where the two patches touch each other

    img_2_cl = snow_ice_fc.reduceToImage(["label"], ee.Reducer.first())  # type: ignore
    bigger_ice = img_2_cl.mask(img_2_cl.select("first").eq(0)).focal_max(2)
    touching_zone = bigger_ice.subtract(img_2_cl.mask(img_2_cl.select("first").gte(9)))
    elev_touch = elevation.addBands(touching_zone, ["first"])
    elevation_snow_line = elev_touch.mask(elev_touch.select("first").lt(-1))
