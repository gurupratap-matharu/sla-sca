import ee

from sla.settings import ING, SRTM


def otsu(hist):
    """
    Return the DN that maximises interclass variance in NIR.
    """

    def get_bss_from_indices(index):
        """
        Compute between sum of squares where each mean partitions the data.
        """

        a_counts = counts.slice(0, 0, index)
        a_count = a_counts.reduce(ee.Reducer.sum(), [0]).get([0])  # type: ignore

        a_means = means.slice(0, 0, index)
        a_mean = (
            a_means.multiply(a_counts)
            .reduce(ee.Reducer.sum(), [0])  # type: ignore
            .get([0])
            .divide(a_count)
        )

        b_count = total.subtract(a_count)
        b_mean = sum.subtract(a_count.multiply(a_mean)).divide(b_count)

        return a_count.multiply(a_mean.subtract(mean).pow(2)).add(
            b_count.multiply(b_mean.subtract(mean).pow(2))
        )

    counts = ee.Array(ee.Dictionary(hist).get("histogram"))  # type: ignore
    means = ee.Array(ee.Dictionary(hist).get("bucketMeans"))  # type: ignore
    size = means.length().get([0])
    total = counts.reduce(ee.Reducer.sum(), [0]).get([0])  # type: ignore
    sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0])  # type: ignore
    mean = sum.divide(total)
    indices = ee.List.sequence(1, size)

    bss = indices.map(get_bss_from_indices)

    # Return the mean value corresponding to the maximum BSS
    return means.sort(bss).get([-1])


def rescale(img, exp, thresholds):
    """
    Applies an expression to an image and linearly rescales it based on threshold value.
    """

    return (
        img.expression(exp, {"img": img})
        .subtract(thresholds[0])
        .divide(thresholds[1] - thresholds[0])
    )


def add_true_hill_shadow(img, elevation, crs_transform):
    """
    TODO
    """

    right_angle = ee.Number(90)
    sun_elevation = ee.Number(img.get("SUN_ELEVATION"))
    sun_azimuth = img.get("SUN_AZIMUTH")
    zenith = right_angle.subtract(sun_elevation)
    shadow_map = (
        ee.Terrain.hillShadow(elevation, sun_azimuth, zenith, 200, True)  # type: ignore
        .focal_min(4)
        .reproject(crs_transform)
    )

    image_4 = img.addBands(shadow_map)

    return image_4


def add_dummy_hill_shadow(img, geometry):
    """
    Add a dummy layer with `-1` value for each pixel.
    """

    shadow_map = ee.Image(-1).clip(geometry).rename(["shadow"])
    image_4 = img.addBands(shadow_map)

    return image_4


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

    # Cast Geometry
    geometry_raw = geometry.geometry()

    # Buffer geometry within 2500m to clip the DSM
    geometry_buffered = geometry_raw.buffer(2500, 5)

    # DEM selection
    dem_glacier = SRTM.select("elevation").rename("AVE_DSM").clip(geometry)

    elevation = dem_glacier.clip(geometry_buffered)

    # Get the CRS information from the original image (Landsat or Sentinel)

    crs_transform = image.select("blue").projection()
    ndsi = image.normalizedDifference(["green", "swir1"]).select("nd").rename("ndsi")
    ndwi = image.normalizedDifference(["green", "nir"]).select("nd").rename("ndwi")

    image_3 = image.addBands(ndsi).addBands(ndwi)

    true_hill_shadow = add_true_hill_shadow(
        img=image_3, elevation=elevation, crs_transform=crs_transform
    )
    dummy_hill_shadow = add_dummy_hill_shadow(img=image_3, geometry=geometry)

    image_with_hill_shadow = ee.Algorithms.If(
        hsboolean, true_hill_shadow, dummy_hill_shadow
    )

    image_4 = ee.Image(image_with_hill_shadow)
    system_time = ee.Number(image.get("system:time_start"))
    hill_shadow_mask = image_4.select("shadow").eq(0).clip(geometry)

    # Start Classification
    blueprint = ee.Image(1).select("constant").rename("classification")
    start = image_4.addBands(blueprint)

    # Get all pixels that are marked as cloudshadows
    cloud_shadow_mask = image.select("cloudShadow").eq(1)
    cloud_hill_shadow_mask = cloud_shadow_mask.add(hill_shadow_mask).gte(1)
    shadows = start.mask(cloud_hill_shadow_mask)

    # Create Shadow Score
    score = ee.Image(1.0)
    score = score.min(
        rescale(img=shadows, exp="img.red + img.green + img.blue", thresholds=[0.5, 2])
    )

    snow_shadow_mask = score.select("constant").gt(0.8)
    unknown_shadow_mask = score.select("constant").lt(0.8)

    snow_shadow = (
        shadows.mask(snow_shadow_mask)
        .select("classification")
        .multiply(6)
        .unmask(-1)
        .clip(geometry)
    )
    unknown_shadow = (
        shadows.mask(unknown_shadow_mask)
        .select("classification")
        .multiply(8)
        .unmask(-1)
        .clip(geometry)
    )

    # Prepare an image without shaded area (shadow is masked)

    image_5 = start.mask(cloud_hill_shadow_mask.Not())

    # Create water mask with ndwi values higher than threshold
    water_mask = image_5.select("ndwi").gt(0.6)

    # Create mask to carry on classification with everything except water
    not_water_mask = image_5.select("ndwi").lt(0.6)

    water_class = start.mask(water_mask)

    shadow = (
        start.mask(water_class.select("blue").gt(0.1))
        .select("classification")
        .multiply(8)
        .unmask(-1)  # Mask everything else with minus 1 (not a class)
    )

    water = (
        start.mask(water_class.select("blue").lt(0.1))
        .select("classification")
        .multiply(2)
        .unmask(-1)  # Mask everything else with minus 1 (not a class)
    )

    image_5 = image_4.mask(not_water_mask)

    snow_ice = image_5.select("ndsi").gt(0.4)
    cloud_debris = image_5.select("ndsi").lt(0.4)

    image_6 = image_5.mask(snow_ice)

    # OTSU Threshold for NIR band

    # compute the histogram of the NIR band

    otsu_image = image_6.clip(geometry).select("nir").multiply(10000)

    histogram = otsu_image.reduceRegion(
        **{
            "reducer": ee.Reducer.histogram(50, 2).combine("mean", None, True),  # type: ignore
            "geometry": geometry,
            "scale": 30,
            "bestEffort": True,
        }
    )

    hist = histogram.get("nir_histogram")
    nir_threshold = otsu(hist=hist).divide(10000)

    # Check whether the nir_threshold is too low or too high
    low, high = nir_threshold.gt(0.41), nir_threshold.lt(0.54)

    is_otsu_high = low.add(high).gt(1)

    fixed_threshold = 0.47

    threshold = ee.Number(
        ee.Algorithms.If(is_otsu_high, nir_threshold, fixed_threshold)
    )

    # Mask Snow and Ice
    snow_mask = image_6.select("nir").gt(threshold)
    snow = start.mask(snow_mask).select("classification").multiply(1).unmask(-1)

    ice_snow_diff = start.mask(image_6.select("nir").lt(threshold))

    score = ee.Image(1.0)

    score = score.min(
        rescale(ice_snow_diff, "img.red + img.green + img.blue", [0.5, 2])
    )

    ice = (
        start.mask(score.select("constant").lt(0.7))
        .select("classification")
        .multiply(0)
        .unmask(1)
    )

    # Mask Cloud and Debris cover

    image_7 = image_5.mask(cloud_debris)
    cloud = (
        start.mask(image_7.select("red").gt(0.6))
        .select("classification")
        .multiply(4)
        .unmask(-1)
    )

    debris = (
        start.mask(image_7.select("red").lt(0.3))
        .select("classification")
        .multiply(3)
        .unmask(-1)
    )

    # Reprojection of all Masks

    water_reprojected = water.reproject(crs_transform)
    cloud_reprojected = cloud.reproject(crs_transform)
    debris_reprojected = debris.reproject(crs_transform)
    snow_reprojected = snow.reproject(crs_transform)
    ice_reprojected = ice.reproject(crs_transform)

    shadow_reprojected = shadow.reproject(crs_transform)
    snow_shadow_reprojected = snow_shadow.reproject(crs_transform)
    unknown_shadow_reprojected = unknown_shadow.reproject(crs_transform)

    # Create one single image with all masks

    classify = (
        water_reprojected.max(cloud_reprojected)
        .max(debris_reprojected)
        .max(snow_reprojected)
        .max(ice_reprojected)
        .max(shadow_reprojected)
        .max(snow_shadow_reprojected)
        .max(unknown_shadow_reprojected)
    ).clip(geometry)

    classified = (
        classify.set("system_time_start", system_time)
        .set("ING_ID", ing_id)
        .set("otsu", threshold)
        .set("deminfo", dem_info)
        .clip(geometry)
        .mask(image.select("blue").gte(0))
    )

    return classified
