import logging
import time

import ee

ee.Initialize()


from sla.classifier import decision_tree  # noqa
from sla.cloud_shadow_mask import add_cloud_shadow  # noqa
from sla.hill_shadow import add_hill_shadow  # noqa
from sla.main_patches import extract_sla_patch  # noqa
from sla.prep import PreProcessor  # noqa

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("main")


def calculate_sla(ing_id="G718255O411666S"):
    """
    Main method which triggers the entire workflow to process a glacier data.
    """

    start = time.time()

    logger.info("starting...")

    start_year = ee.Number(2022)
    end_year = ee.Number(2022)
    doy_start = ee.Number(50)
    doy_end = ee.Number(200)
    cloudiness = ee.Number(50)
    coverage = 50
    hsboolean = 0
    hs_start = ee.Number(50)
    hs_end = ee.Number(200)
    dem = "SRTM"

    logger.info("analysing glacier: %s", ing_id)

    preprocessor = PreProcessor(
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

    logger.info("preprocessing...")
    preprocessed = preprocessor.execute()

    logger.info("processing cloud shadow mask...")
    cloud_shadow_collection = preprocessed.map(add_cloud_shadow)

    logger.info("processing hill shadow...")
    hill_shadow_collection = add_hill_shadow(
        image_collection=cloud_shadow_collection, hs_start=hs_start, hs_end=hs_end
    )

    logger.info("classifying with decision tree...")
    classified_collection = hill_shadow_collection.map(decision_tree)

    logger.info("extracting snow line altitude (sla)...")
    sla = classified_collection.map(extract_sla_patch)

    end = time.time()
    logger.info("All Done...💅🏻💫💖!")
    logger.info(f"It took {end-start:.2f} seconds")

    return sla


if __name__ == "__main__":
    ING_ID = "G718255O411666S"  # <-- ONLY change this to try a different glacier
    SLA = calculate_sla(ing_id=ING_ID)
    # write_to_local(response=SLA, filename="dump/sla.json")
