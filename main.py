import json
import logging

import ee

ee.Initialize()

from cloud_shadow_mask import add_cloud_shadow  # noqa
from prep import PreProcessor  # noqa

logger = logging.getLogger(__name__)


def write_to_local(response, filename):
    """
    Writes a GEE object to the local filesystem in JSON format.

    Be cautious as this method will call the getInfo() method on the GEE object to
    retrieve all the results from the EE server. So it can be computationally expensive
    and block all further execution.
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
preprocessed = preprocessed.filterMetadata(
    "system:time_start", "not_equals", 1443176819706
)
print("processing cloud shadow mask...")
preprocessed_with_cloud_shadows = preprocessed.map(add_cloud_shadow)


write_to_local(
    response=preprocessed_with_cloud_shadows, filename="dump/cloud_shadow.json"
)

print("All Done...!")
