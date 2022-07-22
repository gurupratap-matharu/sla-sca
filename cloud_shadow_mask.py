""""
Script to identify pixels in image which are under a cloud shadow area.
We mark these areas for correct classification and detection of glacier area

The script exposes a CloudShadowMask class which holds various utility helper methods defined within.
This `CloudShadowMask` object is supposed to be used in a Map Function over the entire ImageCollection
to mask relevant pixels.
"""


import logging
import math

import ee

logger = logging.getLogger(__name__)

FC_PATH = "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve"
FC = ee.FeatureCollection(FC_PATH)


def rescale(img, exp, thresholds):
    return (
        img.expression(exp, {"img": img})
        .subtract(thresholds[0])
        .divide(thresholds[1] - thresholds[0])
    )


def add_cloud_score(img):
    ndmi = img.normalizedDifference(["nir", "swir1"])
    ndsi = img.normalizedDifference(["green", "swir1"])

    score = ee.Image(1)
    score = score.min(rescale(img, "img.blue", [0.1, 0.5]))
    score = score.min(rescale(img, "img.red + img.green + img.blue", [0.2, 0.8]))
    score = score.min(rescale(ndmi, "img", [-0.1, 0.1]))
    score = score.min(rescale(ndsi, "img", [0.4, 0.1]))
    score = score.multiply(100).byte()

    return img.addBands(score.rename("cloudScore"))


def project_shadows(
    image, cloud_mask, cloud_heights, contract_pixels, dilate_pixels, infrared_threshold
):
    def find_shadow(img):
        # TODO
        return img

    right_angle = ee.Number(90)
    sun_elevation = image.get("SUN_ELEVATION")
    mean_azimuth = image.get("SUN_AZIMUTH")
    mean_zenith = right_angle.subtract(sun_elevation)

    # Find dark pixels in the image
    dark_pixels = (
        image.select(["nir", "swir1", "swir2"])
        .reduce(ee.Reducer.sum())
        .lt(infrared_threshold)
    )

    # Get Scale of image
    nominal_scale = cloud_mask.projection().nominalScale()

    # Find where cloud shadows should be based on Solar Geometry and convert it to radians
    azimuth_radians = (
        ee.Number(mean_azimuth)
        .multiply(math.pi)
        .divide(180.0)
        .add(ee.Number(0.5).multiply(math.pi))
    )
    zenith_radians = (
        ee.Number(0.5)
        .multiply(math.pi)
        .subtract(ee.Number(mean_zenith).multiply(math.pi).divide(180.0))
    )

    # Find the shadows
    shadows = cloud_heights.map(find_shadow)
    shadow_mask = ee.ImageCollection.fromImages(shadows).max()
    # shadow_mask = shadow_mask.and(cloud_mask.not())
    # shadow_mask = shadow_mask.focal_min(contract_pixels).focal_max(dilate_pixels)
    # shadow_mask = shadow_mask.and(dark_pixels)

    image = image.addBands(shadow_mask.rename(["cloudShadow"]))

    return image


def add_cloud_shadow(image):
    ing_id = image.get("ING_ID")
    cloud_threshold = 20
    cloud_heights = ee.List.sequence(0, 500, 25)
    contract_pixels = 10
    dilate_pixels = 15
    infrared_threshold = 0.45
    geometry = FC.filterMetadata("ID_local", "equals", ing_id)

    # Get cloud score
    image = add_cloud_score(image)
    cloud_score = image.select("cloudScore")

    # Mask the cloud score
    cloud_mask = (
        cloud_score.gt(cloud_threshold)
        .focal_min(contract_pixels)
        .focal_max(dilate_pixels)
    )
    # CloudMasked  = image.updateMask(cloud_mask.not()) # TODO check this.... is this needed?

    # Projecting the cloud shadow
    cloud_shadow_masked = project_shadows(
        image,
        cloud_mask,
        cloud_heights,
        contract_pixels,
        dilate_pixels,
        infrared_threshold,
    )
    return cloud_shadow_masked.clip(geometry)


if __name__ == "__main__":
    pass
