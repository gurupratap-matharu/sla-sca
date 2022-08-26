import json
import logging
import time

import ee

ee.Initialize()

from sla.classifier import decision_tree  # noqa
from sla.cloud_shadow_mask import add_cloud_shadow  # noqa
from sla.hill_shadow import add_hill_shadow  # noqa
from sla.prep import PreProcessor  # noqa

logger = logging.getLogger(__name__)

start = time.time()

print("starting...")


def write_to_local(response, filename):
    """
    Writes a GEE object to the local filesystem in JSON format.

    Be cautious as this method will call the getInfo() method on the GEE
    object to retrieve all the results from the EE server.
    So it can be computationally expensive and block all further execution.
    """

    with open(filename, "w", encoding="utf-8") as f:
        f.write(json.dumps(response.getInfo()))

    print(f"Results written to {filename}")


ing_id = "G718255O411666S"
start_year = ee.Number(2022)
end_year = ee.Number(2022)
doy_start = ee.Number(50)
doy_end = ee.Number(200)
cloudiness = ee.Number(50)
coverage = 50
hsboolean = 0
dem = "SRTM"


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

print("preprocessing...")
preprocessed = preprocessor.execute()

print("processing cloud shadow mask...")
preprocessed = preprocessed.map(add_cloud_shadow)

print("processing hill shadow...")
preprocessed = add_hill_shadow(image_collection=preprocessed)


print("initiating classifier...")
map_collection = preprocessed.map(decision_tree)

print("writing to local...")
write_to_local(response=map_collection, filename="dump/classified.json")

end = time.time()
print("All Done...!")
print(f"It took {end-start:.2f} seconds")
