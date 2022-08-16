import imghdr

import ee

ALOS = ee.Image("JAXA/ALOS/AW3D30/V2_2")
SRTM = ee.Image("USGS/SRTMGL1_003")
ING = ee.FeatureCollection(
    "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve"
)


def hill_shadow_yes():
    pass


def hill_shadow_no():
    pass


def decision_tree(image):
    """
    Classifies each image... TODO
    """
    # Get input values from image to be classified

    # Check if hill shadow boolean added to each image
    hsboolean = image.get("hsboolean")

    # Get GLIMS ID
    ing_id = image.get("ING_ID")
    dem_info = image.get("deminfo")

    # Get shape file
    geometry = ING.filterMetadata("ID_local", "equals", ing_id)

    geometry_raw = geometry.geometry()
    geometry_buffered = geometry_raw.buffer(2500, 5)

    # DEM selection
    dem_1 = SRTM.select("elevation").rename("AVE_DSM").clip(geometry)

    dem_2 = ALOS.select("AVE_DSM").clip(geometry)

    # Check if dem selection == 'ALOS'
    dem_selector = ee.Algorithms.IsEqual(ee.String(dem_info), ee.String("ALOS"))

    # Conditionally select the DEM to be used
    # TODO CHECK this.. We should avoid If statements!
    dem_glacier = ee.Image(ee.Algorithms.If(dem_selector, dem_2, dem_1))

    elevation = dem_glacier.clip(geometry_buffered)

    # Get the CRS information from the original image (Landsat or Sentinel)

    crs_transform = image.select("blue").projection()
    ndsi = image.normalizedDifference(["green", "swir1"]).select("nd").rename("ndsi")
    ndwi = image.normalizedDifference(["green", "nir"]).select("nd").rename("ndwi")

    image2 = image.addBands(ndsi)
    image3 = image2.addBands(ndwi)

    add_hill_shadow = ee.Algorithms.If(
        hsboolean, hill_shadow_yes(), hill_shadow_no()
    )  # check if these functions should be called here or only referred

    image4 = ee.Image(add_hill_shadow)
    system_time = ee.Number(image.get("system:time_start"))
    hill_shadow_mask = image4.select("shadow").eq(0).clip(geometry)

    return image
