import ee

ING = ee.FeatureCollection(
    "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve"
)
ALOS = ee.Image("JAXA/ALOS/AW3D30/V2_2")
SRTM = ee.Image("USGS/SRTMGL1_003")
