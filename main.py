import json

import ee

ee.Initialize()
import logging

from cloud_shadow_mask import add_cloud_shadow
from prep import PreProcessor

logger = logging.getLogger(__name__)


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
print("processing cloud shadowj mask...")
preprocessed_with_cloud_shadows = preprocessed.map(add_cloud_shadow)


# Write to local file
response = preprocessed_with_cloud_shadows.getInfo()
FILENAME = "cloud_shadow_masked.json"

with open(FILENAME, "w", encoding="utf-8") as f:
    f.write(json.dumps(response))

    print(f"Results written to {FILENAME}")

print("All Done...!")
