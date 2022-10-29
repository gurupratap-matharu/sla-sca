import logging
import time

import ee

ee.Initialize()


from sla.classifier import decision_tree  # noqa
from sla.cloud_shadow_mask import add_cloud_shadow  # noqa
from sla.hill_shadow import add_hill_shadow  # noqa
from sla.main_patches import extract_sla_patch  # noqa
from sla.prep import PreProcessor  # noqa
from sla.utils import write_to_local

logging.basicConfig(level=logging.INFO)

logger = logging.getLogger("main")


def main():
    """
    Main method which triggers the entire workflow to process a glacier data.
    """

    start = time.time()

    logger.info("starting...")

    ING_IDS = ("G718255O411666S",)
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

    for ing_id in ING_IDS:
        logger.info("analysing id: %s", ing_id)

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
        preprocessed = preprocessed.map(add_cloud_shadow)

        logger.info("processing hill shadow...")
        preprocessed = add_hill_shadow(
            image_collection=preprocessed, hs_start=hs_start, hs_end=hs_end
        )

        logger.info("initiating classifier...")
        map_collection = preprocessed.map(decision_tree)

        logger.info("initiating SLA extraction...")
        map_collection = map_collection.map(extract_sla_patch)

        # logger.info("writing to local...")
        # write_to_local(response=map_collection, filename="dump/final.json")

    end = time.time()
    logger.info("All Done...!")
    logger.info(f"It took {end-start:.2f} seconds")


if __name__ == "__main__":
    main()
