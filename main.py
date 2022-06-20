import ee

ee.Initialize()
from prep import PreProcessor

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

    res = obj.execute()
    print(res)
