"""
Preparation Module for SLA - SCA

This script composes of various helper functions but its main role is to prepare and provide an image collection.
The script exposes a high level function which does all the heavy lifting.

The composite image is generated with this script and used in main code for further processing.
"""
import json

import ee

from sla.settings import ING, L5, L7, L8, S2

ee.Initialize()


FILENAME = "dump/prep.json"


class PreProcessor:
    """
    Main Preprocess class for all Sensor Data

    Filters data from all Sensors (L5, L7, L8 and S2) via date and other parameters
    and builds a single ImageCollection needed for main processing.
    """

    def __init__(
        self,
        ing_id,
        start_year,
        end_year,
        cloudiness,
        coverage,
        doy_start,
        doy_end,
        hsboolean,
        dem,
    ):
        self.ing_id = ing_id

        self.start_date = ee.Date.fromYMD(start_year, 1, 1)
        self.end_date = ee.Date.fromYMD(end_year, 12, 31)
        self.calendar_filter = ee.Filter.calendarRange(
            doy_start, doy_end, "day_of_year"
        )

        self.cloud_filter = ee.Filter.gt("noclouds", cloudiness)
        self.coverage_filter = ee.Filter.gt("arearatio", coverage)

        self.hsboolean = hsboolean
        self.dem = dem

        self.l5_bands = ee.List(["B1", "B2", "B3", "B4", "B5", "B6", "B7", "BQA"])
        self.l5_band_names = ee.List(
            ["blue", "green", "red", "nir", "swir1", "tir1", "swir2", "BQA"]
        )

        self.l7_bands = ee.List(
            ["B1", "B2", "B3", "B4", "B5", "B6_VCID_1", "B7", "BQA"]
        )

        self.l7_band_names = ee.List(
            ["blue", "green", "red", "nir", "swir1", "tir1", "swir2", "BQA"]
        )

        self.l8_bands = ee.List(
            ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B10", "B11", "BQA"]
        )

        self.l8_band_names = ee.List(
            [
                "cb",
                "blue",
                "green",
                "red",
                "nir",
                "swir1",
                "swir2",
                "tir1",
                "tir2",
                "BQA",
            ]
        )

        self.s2_bands = ee.List(
            [
                "B1",
                "B2",
                "B3",
                "B4",
                "B5",
                "B6",
                "B7",
                "B8",
                "B8A",
                "B9",
                "B10",
                "B11",
                "B12",
                "QA60",
            ]
        )

        self.s2_band_names = ee.List(
            [
                "cb",
                "blue",
                "green",
                "red",
                "re1",
                "re2",
                "re3",
                "nir",
                "re4",
                "vapor",
                "cirrus",
                "swir1",
                "swir2",
                "BQA",
            ]
        )

        self.geometry = None
        self.geometry_raw = None
        self.geometry_buffered = None
        self.area_glacier_IMG = None

        self.l5_filtered = None
        self.l7_filtered = None
        self.l8_filtered = None
        self.s2_filtered = None

    def execute(self):
        """
        Main method to start preprocessing for all sensor data.
        """

        self.geometry = ING.filterMetadata("ID_local", "equals", self.ing_id)
        self.geometry_raw = self.geometry.geometry()
        self.geometry_buffered = self.geometry_raw.buffer(2500, 5)

        glacier_area = (
            ee.Image(1).clip(self.geometry).select("constant").rename("green")
        )

        # TODO check if areatotal is needed?
        # areatotal = self.geometry.first().get("Area")
        self.area_glacier_IMG = self.area_calc_rast(glacier_area)

        self.l5_filtered = self.filter_l5()
        self.l7_filtered = self.filter_l7()
        self.l8_filtered = self.filter_l8()
        self.s2_filtered = self.filter_s2()

        collection = (
            self.l5_filtered.merge(self.l7_filtered)
            .merge(self.l8_filtered)
            .merge(self.s2_filtered)
        )

        collection = collection.map(self.remove_duplicates)

        collection = collection.filterMetadata(
            "system:time_start", "not_equals", 1443176819706
        )

        return collection

    def filter_l5(self):
        """
        Filter L5 sensor data based on date and geometry and applies various filters
        like cloud and coverage.
        """

        return (
            L5.filterBounds(self.geometry)
            .filterDate(self.start_date, self.end_date)
            .filter(self.calendar_filter)
            .map(lambda image: image.select(self.l5_bands).rename(self.l5_band_names))
            .map(self.add_quality_info)
            .filter(self.cloud_filter)
            .filter(self.coverage_filter)
            .map(self.add_albedo)
            .sort("system:time_start")
        )

    def filter_l7(self):
        """

        Filter L7 sensor data based on date and geometry and applies various filters
        like cloud and coverage.
        """

        return (
            L7.filterBounds(self.geometry)
            .filterDate(self.start_date, self.end_date)
            .filter(self.calendar_filter)
            .map(lambda image: image.select(self.l7_bands).rename(self.l7_band_names))
            .map(self.add_quality_info)
            .filter(self.cloud_filter)
            .filter(self.coverage_filter)
            .map(self.add_albedo)
            .sort("system:time_start")
        )

    def filter_l8(self):
        """

        Filter L5 sensor data based on date and geometry and applies various filters
        like cloud and coverage.
        """

        return (
            L8.filterBounds(self.geometry)
            .filterDate(self.start_date, self.end_date)
            .filter(self.calendar_filter)
            .map(lambda image: image.select(self.l8_bands).rename(self.l8_band_names))
            .map(self.add_quality_info)
            .filter(self.cloud_filter)
            .filter(self.coverage_filter)
            .map(self.add_albedo)
            .sort("system:time_start")
        )

    def filter_s2(self):
        """
        Filter S2 sensor data based on date and geometry and applies various filters
        like cloud and coverage.
        """

        return (
            S2.filterBounds(self.geometry)
            .filterDate(self.start_date, self.end_date)
            .filter(self.calendar_filter)
            .map(lambda image: image.select(self.s2_bands).rename(self.s2_band_names))
            .map(self.sentinel_add_quality_info)
            .filter(self.cloud_filter)
            .filter(self.coverage_filter)
            .map(self.scale_s2_pixels)
            .map(self.add_albedo)
            .sort("system:time_start")
        )

    def rescale(self, img, exp, thresholds):
        """
        Applies an expression to an image and linearly rescales it based on threshold value.
        """

        return (
            img.expression(exp, {"img": img})
            .subtract(thresholds[0])
            .divide(thresholds[1] - thresholds[0])
        )

    def cloud_score(self, img):
        """
        Compute the cloud score for masking. It expects the input image to have common band names.

        Cloud Score is only for the nocloud-ratio and has no influence on the classification.
        If the cloud-ratio detected with this score too high then it's possible that the images are filtered
        within the next step but all images are included in the collection.
        """

        # Compute several indicators of cloudiness and take the minimum of them.
        score = ee.Image(1.0)

        # Clouds are reasonably bright in the blue band.
        score = score.min(self.rescale(img, "img.blue", [0.1, 0.3]))

        # Clouds are reasonably bright in all visible bands.
        score = score.min(
            self.rescale(img, "img.red + img.green + img.blue", [0.2, 0.8])
        )

        # Clouds are reasonably bright in all infrared bands.
        score = score.min(
            self.rescale(img, "img.nir + img.swir1 + img.swir2", [0.3, 0.8])
        )

        # However, clouds are not snow.
        ndsi = img.normalizedDifference(["green", "swir1"])
        score2 = score.min(self.rescale(ndsi, "img", [0.7, 0.6]))
        score3 = ee.Image(1).subtract(score2).select([0], ["cloudscore"])

        # add band with cloud score per pixel and select cloud pixels to mask original image
        img1_1 = img.addBands(score3)
        img1_2 = img1_1.select("cloudscore").gt(0.7)
        img1_3 = img.mask(img1_2)

        return img1_3

    def area_calc_rast(self, img):
        """
        Calculate the area of an image (area of a band without masked pixels)
        """

        # select not masked pixels
        img1 = img.select("green").neq(0)
        fc1 = img1.reduceToVectors(
            **{
                "reducer": ee.Reducer.countEvery(),  # type: ignore
                "geometry": self.geometry,
                "scale": 30,
                "maxPixels": 1e10,
                "bestEffort": True,
                "crs": "EPSG:4326",
            }
        )

        area = fc1.reduceColumns(
            **{
                "reducer": ee.Reducer.sum(),  # type: ignore
                "selectors": ["count"],
            }
        )

        return ee.Number(area.get("sum")).divide(1000)

    def area_calc(self, image):
        """
        Calculates the green colored area in an image
        """

        prep = image.select("green")
        pixelarea = prep.multiply(ee.Image.pixelArea())

        return pixelarea.reduceRegion(
            **{
                "reducer": ee.Reducer.sum(),  # type: ignore
                "geometry": self.geometry,
                "scale": 30,
                "bestEffort": True,
                "maxPixels": 1e10,
            }
        ).getNumber("green")

    def add_albedo(self, img):
        """
        Add Albedo to each image based on LIANG (2000) - Narrowband to broadband conversions
        of land surface albedo: I Algorithms
        """

        albedo = (
            img.select("blue")
            .multiply(0.356)
            .add(img.select("red").multiply(0.130))
            .add(img.select("nir").multiply(0.373))
            .add(img.select("swir1").multiply(0.085))
            .add(img.select("swir2").multiply(0.072))
            .subtract(0.0018)
            .rename(["albedo"])
        ).divide(1.016)

        return img.addBands(albedo)

    def cloud_nocloud_ratio(self, data):
        """
        Calculate the cloud/no-cloud ratio per image.
        """

        clipped_data = data.clip(self.geometry)
        cloudfreefunc = self.cloud_score(clipped_data)

        # calculate the area for the image with clouds and when the detected clouds are masked
        countcf, countnorm = self.area_calc(cloudfreefunc), self.area_calc(clipped_data)

        return countcf.divide(countnorm).multiply(100)

    def add_quality_info(self, data):
        """
        Add all the 'Quality-Information' to the image
        """

        # TODO: VEER CHECK THIS FUNCTION. DO WE NEED RGI
        # TACKLE GLOBAL VARIABLES

        clip = data.clip(self.geometry_buffered)
        quality = self.cloud_nocloud_ratio(clip)
        image = clip.set("noclouds", quality)
        sens = clip.get("SPACECRAFT_ID")

        # calculate the area for the image to detect, if there are parts without information
        area_rgi_outline = self.area_calc_rast(clip)  # TODO DO WE NEED THIS? ASK
        # RGI:
        # ING: inventario nacional de glaciares de Argentina

        # calculate ratio of coverage from one image over the glacier
        arearatio = area_rgi_outline.divide(self.area_glacier_IMG).multiply(100)

        return (
            image.set("ING_ID", self.ing_id)
            .set("areaRGIoutline", area_rgi_outline)
            .set("areaIMGglacier", self.area_glacier_IMG)
            .set("arearatio", arearatio)
            .set("SENSOR", sens)
            .set("hsboolean", self.hsboolean)
            .set("deminfo", self.dem)
        )

    def cloud_filter_sentinel(self, data):
        """
        Normalizes the cloud regions in sentinel images data
        """

        clipped_data = data.clip(self.geometry)

        countnorm = self.area_calc(clipped_data)
        cloudfreefunc = self.cloud_score(clipped_data)
        countcf = self.area_calc(cloudfreefunc)

        return countcf.divide(countnorm).multiply(100)

    def sentinel_add_quality_info(self, data):
        """
        Add quality information for each sentinel image.
        """

        clip = data.clip(self.geometry_buffered)
        quality = self.cloud_filter_sentinel(clip)
        area_rgi_outline = self.area_calc_rast(clip)
        area_ing_outline = self.area_calc_rast(clip)
        arearatio = area_rgi_outline.divide(self.area_glacier_IMG).multiply(100)
        rightangle = ee.Number(90)
        sun_az = ee.Number(clip.get("MEAN_SOLAR_AZIMUTH_ANGLE"))
        sun_el_zenith = ee.Number(clip.get("MEAN_SOLAR_ZENITH_ANGLE"))
        sens = clip.get("SPACECRAFT_NAME")
        sun_el = rightangle.subtract(ee.Number(sun_el_zenith))
        image = clip.set("noclouds", quality)

        return (
            image.set("ING_ID", self.ing_id)
            .set("areaINGoutline", area_ing_outline)
            .set("areaIMGglacier", self.area_glacier_IMG)
            .set("arearatio", arearatio)
            .set("SUN_AZIMUTH", sun_az)
            .set("SUN_ELEVATION", sun_el)
            .set("SENSOR", sens)
            .set("hsboolean", self.hsboolean)
            .set("deminfo", self.dem)
        )

    def scale_s2_pixels(self, image):
        """
        Scales the images from Sentinel sensor to a ratio between 0 and 1
        """

        custom_bands = [
            "cb",
            "blue",
            "green",
            "red",
            "re1",
            "re2",
            "re3",
            "nir",
            "re4",
            "vapor",
            "cirrus",
            "swir1",
            "swir2",
        ]

        rescaled_image = image.select(custom_bands).divide(10000)
        return image.addBands(srcImg=rescaled_image, overwrite=True)

    def remove_duplicates(self, image):
        """
        Filters out all duplicate images in an image collection.
        """

        # TODO: FIND OUT LOGIC FOR THIS. ORIGINAL IMPLEMENTATION SEEMS AMBIGIOUS
        return image


if __name__ == "__main__":

    ing_id = "G718255O411666S"
    start_year = ee.Number(2022)
    end_year = ee.Number(2022)
    doy_start = ee.Number(50)
    doy_end = ee.Number(200)
    cloudiness = ee.Number(50)
    coverage = 50
    hsboolean = 0
    dem = "SRTM"

    obj = PreProcessor(
        ing_id=ing_id,
        start_year=start_year,
        end_year=end_year,
        doy_start=doy_start,
        doy_end=doy_end,
        cloudiness=cloudiness,
        coverage=coverage,
        hsboolean=hsboolean,
        dem=dem,
    )

    print("executing...")
    response = obj.execute()
    print("getting info...")
    response = response.getInfo()

    with open(FILENAME, "w", encoding="utf-8") as f:
        f.write(json.dumps(response))

    print(f"Results written to {FILENAME}")
    print("All done...!")
