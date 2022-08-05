import ee

HS_START = ee.Number(90)
HS_END = ee.Number(100)


def add_hill_shadow(image_collection, hs_start=HS_START, hs_end=HS_END):

    # Split the collection into parts to write a 0 or 1
    hsy_filter = ee.Filter.calendarRange(hs_start, hs_end, "day_of_year")
    hsn_filter = ee.Filter.calendarRange(hs_end, hs_start, "day_of_year")

    # Add a Yes to the 'hsboolean' property (hill shadow boolean)
    img_coll_hsy = image_collection.filter(hsy_filter).map(
        lambda img: img.set("hsboolean", 1)
    )

    # Add a No to the 'hsboolean' property (hill shadow boolean)
    img_coll_hsn = image_collection.filter(hsn_filter).map(
        lambda img: img.set("hsboolean", 0)
    )

    # Merge the two parts and sort them in ascending order

    prep_final = ee.ImageCollection(
        img_coll_hsn.merge(img_coll_hsy).sort("system:time_start")
    )

    return prep_final
