"""
Preparation Module for SLA - SCA

This script composes of various helper functions but its main role is to prepare and provide an image collection. The script exposes a high level
function which does all the heavy lifting.

The composite image is generated with this script and used in main code for further processing.
"""

import ee

ING = ee.FeatureCollection(
    "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve")
L5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA")
L7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_TOA")
L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
S2 = ee.ImageCollection("COPERNICUS/S2")

"""
Test data....

ing_id = 'G718255O411666S'
startyear =           ee.Number(2022) #
endyear =             ee.Number(2022) # Año de fin del estudio
filterDOYstart =      ee.Number(50)   # Start of Range: use DOY-Number to filter the range
filterDOYend =        ee.Number(200)  # End of Range: DOY-Number to filter the range
cloudiness =          ee.Number(50)   #% (Area over the Glacier that is not covered with Clouds)
coverage =            50             #percentage of the glacier area, that has to be cove


https://code.earthengine.google.com/?asset=users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve
"""


def imageCollection(ing_id, start_year, end_year, cloudiness, coverage, doy_start, doy_end, hsboolean, dem):

    l5_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'BQA'])
    l5_band_names = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'tir1', 'swir2', 'BQA'])

    l7_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_1', 'B7', 'BQA'])
    l7_band_names = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'tir1', 'swir2', 'BQA'])

    l8_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'B11', 'BQA'])
    l8_band_names = ee.List(['cb', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir1', 'tir2', 'BQA'])

    s2_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12', 'QA60'])
    s2_band_names = ee.List(['cb', 'blue', 'green', 'red', 're1', 're2', 're3',
                            'nir', 're4', 'vapor', 'cirrus', 'swir1', 'swir2', 'BQA'])

    start_date = ee.Date.fromYMD(start_year, 1, 1)
    end_date = ee.Date.fromYMD(end_year, 12, 31)
    calendar_filter = ee.Filter.calendarRange(doy_start, doy_end, 'day_of_year')

    cloud_filter = ee.Filter.gt('noclouds', cloudiness)
    coverage_filter = ee.Filter.gt('arearatio', coverage)

    geometry = ING.filterMetadata('ID_local', 'equals', ing_id)
    geometry_raw = geometry.geometry()
    geometry_buffered = geometry_raw.buffer(2500, 5)

    # Calculate glacier area
    areatotal1 = ee.Image(1).clip(geometry).select('constant').rename('green')

    areatotal = geometry.first().get('Area')
    areaIMGglacier1 = area_calc_rast(areatotal1)

    L5_FILTERED = L5.filterBounds(geometry) \
                    .filterDate(start_date, end_date) \
                    .filter(calendar_filter)

    L7_FILTERED = L7.filterBounds(geometry) \
                    .filterDate(start_date, end_date) \
                    .filter(calendar_filter)

    L8_FILTERED = L8.filterBounds(geometry) \
                    .filterDate(start_date, end_date) \
                    .filter(calendar_filter)


def rescale(img, exp, thresholds):
    """
    Applies an expression to an image and linearly rescales it based on threshold value.
    """

    return img.expression(exp, {'img': img}).subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])


def cloud_score(img):
    """
    Compute the cloud score for masking. It expects the input image to have common band names.

    Cloud Score is only for the nocloud-ratio and has no influence on the classification.
    If the cloud-ratio detected with this score too high then it's possible that the images are filtered
    within the next step but all images are included in the collection.
    """

    # Compute several indicators of cloudyness and take the minimum of them.
    score = ee.Image(1.0)

    # Clouds are reasonably bright in the blue band.
    score = score.min(rescale(img, 'img.blue', [0.1, 0.3]))

    # Clouds are reasonably bright in all visible bands.
    score = score.min(rescale(img, 'img.red + img.green + img.blue', [0.2, 0.8]))

    # Clouds are reasonably bright in all infrared bands.
    score = score.min(rescale(img, 'img.nir + img.swir1 + img.swir2', [0.3, 0.8]))

    # However, clouds are not snow.
    ndsi = img.normalizedDifference(['green', 'swir1'])
    score2 = score.min(rescale(ndsi, 'img', [0.7, 0.6]))
    score3 = ee.Image(1).subtract(score2).select([0], ['cloudscore'])

    # add band with cloud score per pixel and select cloud pixels to mask original image
    img1_1 = img.addBands(score3)
    img1_2 = img1_1.select('cloudscore').gt(0.7)
    img1_3 = img.mask(img1_2)

    return img1_3


def area_calc_rast(img):
    """
    Calculate the area of an image (area of a band without masked pixels)
    """

    # select not masked pixels
    img1 = img.select('green').neq(0)
    fc1 = img1.reduceToVectors({
        'reducer': ee.Reducer.countEvery(),
        'geometry': geometry,
        'scale': 30,
        'maxPixels': 1e10,
        'bestEffort': True,
        'crs': 'EPSG:4326'
    })

    area = fc1.reduceColumns({
        'reducer': ee.Reducer.sum(),
        'selectors': ['count'],
    })

    return ee.Number(area.get('sum')).divide(1000)


def area_calc(image):

    prep = image.select('green')
    pixelarea = prep.multiply(ee.Image.pixelArea())

    return pixelarea.reduceRegion({
        'reducer': ee.Reducer.sum(),
        'geometry': geometry,
        'scale': 30,
        'bestEffort': True,
        'maxPixels': 1e10
    }).getNumber('green')


def add_albedo(img):
    """
    Add Albedo to each image based on LIANG (2000) - Narrowband to broadband conversions
    of land surface albedo: I Algorithms
    """

    albedo = (((img.select('blue').multiply(0.356))
                    .add(img.select('red').multiply(0.130))
                    .add(img.select('nir').multiply(0.373))
                    .add(img.select('swir1').multiply(0.085))
                    .add(img.select('swir2').multiply(0.072))
                    .subtract(0.0018).rename(['albedo'])).divide(1.016))

    return img.addBands(albedo)


def cloud_nocloud_ratio(data):
    """
    Calculate the cloud/no-cloud ratio per image.
    """

    input = data.clip(geometry)
    cloudfreefunc = cloud_score(input)

    # calculate the area for the image with clouds and when the detected clouds are masked
    countcf, countnorm = area_calc(cloudfreefunc), area_calc(input)

    return countcf.divide(countnorm).multiply(100)


def add_quality_info(data):
    """
    Add all the 'Quality-Information' to the image
    """

    # TODO: VEER CHECK THIS FUNCTION. DO WE NEED RGI
    # TACKLE GLOBAL VARIABLES

    clip = data.clip(geometrybuffered)
    quality = cloud_nocloud_ratio(clip)
    image = clip.set('noclouds', quality)
    sens = clip.get('SPACECRAFT_ID')

    # calculate the area for the image to detect, if there are parts without information
    areaRGIoutline = area_calc_rast(clip)  # TODO DO WE NEED THIS? ASK
    areaIMGglacier = areaIMGglacier1  # areacalcrast(areatotal1)

    # calculate ratio of coverage from one image over the glacier
    arearatio = areaRGIoutline.divide(areaIMGglacier).multiply(100)
    return (image.set('ING_ID', ing_id)
                    .set('areaRGIoutline', areaRGIoutline)
                    .set('areaIMGglacier', areaIMGglacier)
                    .set('arearatio', arearatio)
                    .set('SENSOR', sens)
                    .set('hsboolean', hsboolean)
                    .set('deminfo', dem))


def cloud_filter_sentinel(data):

    input = data.clip(geometry)
    cloudfreefunc = cloud_score(input)
    countcf = area_calc(cloudfreefunc)
    countnorm = area_calc(input)
    return countcf.divide(countnorm).multiply(100)


def s_add_quality_info(data):

    clip = data.clip(geometrybuffered)
    quality = cloud_filter_sentinel(clip)
    areaRGIoutline = area_calc_rast(clip)
    areaIMGglacier = areaIMGglacier1
    arearatio = areaRGIoutline.divide(areaIMGglacier).multiply(100)
    rightangle = ee.Number(90)
    sun_az = ee.Number(clip.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
    sun_el_zenith = ee.Number(clip.get('MEAN_SOLAR_ZENITH_ANGLE'))
    sens = (clip.get('SPACECRAFT_NAME'))
    sunEL = rightangle.subtract(ee.Number(sun_el_zenith))
    image = clip.set('noclouds', quality)

    return (image.set('ING_ID', ing_id)
        .set('areaINGoutline', areaINGoutline)
        .set('areaIMGglacier', areaIMGglacier)
        .set('arearatio', arearatio)
        .set('SUN_AZIMUTH', sun_az)
        .set('SUN_ELEVATION', sunEL)
        .set('SENSOR', sens)
        .set('hsboolean', hsboolean)
        .set('deminfo', dem))


def func_owr(img)  # create image collection:
                  img1 = img.select(l5_bands).rename(l5_band_names);  # rename all images with the correct names
                  img2 = img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_owr) \
                .map(addquainfo) \
                .filter(cloud_filter) \
                .filter(coverage_filter) \
                .map(addAlbedo)


def func_xbh(img):
                  img1 = img.select(l7_bands).rename(l7_band_names)
                  img2 = img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                  .map(func_xbh) \
                .map(addquainfo) \
                .filter(cloud_filter) \
                .filter(coverage_filter) \
                .map(addAlbedo)


def func_tfu(img):
                  img1 = img.select(l8_bands).rename(l8_band_names)
                  img2 = img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_tfu) \
                .map(addquainfo) \
                .filter(cloud_filter) \
                .filter(coverage_filter) \
                .map(addAlbedo)


sentineldata = s2.filterBounds(geometry) \
                .filterDate(startdate, enddate).filter(ee.Filter.calendarRange(filterDOYstart, filterDOYend, 'day_of_year'))


def func_ngk(img):
                  img1 = img.select(s2_bands).rename(s2_band_names)
                  img2 = img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_ngk) \
                .map(saddquainfo) \
                .filter(cloud_filter) \
                .filter(coverage_filter)


def func_jty(img):
                  # scale the pixelvalues to a ratio between 0 and 1
                  img1 = img \
                    .select((['cb', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'vapor', 'cirrus', 'swir1', 'swir2'])) \
                    .divide(10000)
                  return img.addBands(img1, None, True) \
                .map(func_jty
).map(addAlbedo)


).map(addAlbedo)


###########   FILTER DUBLICATES WITH SAME DATE    #########################
# sort all collection (for each sensor), that the dublicate filter works
landsat5datapre=L5_FILTERED.sort('system:time_start')
landsat7datapre=landsat7data.sort('system:time_start')
landsat8datapre=landsat8data.sort('system:time_start')
sentineldatapre=sentineldata.sort('system:time_start')
# Function to detect dublicates
def startdublicatfiltering(imgcollection):


condition1=imgcollection.size()


def normal():


list=imgcollection.toList(imgcollection.size())
image=ee.Image(list.get(0))
dummyimage=ee.Image(1).set('system:time_start', 0)
# Add in the end of the list a dummy image
list=list.add(dummyimage)


def detect_dublicates(image):
  isdublicate=ee.String("")
  number=list.indexOf(image)
  image1=ee.Image(list.get(number.add(1)))
  # Compare the image(0) in the ImageCollection with the image(1) in the List
  date1=image.date().format("Y-M-d")
  date2=image1.date().format("Y-M-d")
  cond=ee.Algorithms.IsEqual(date1, date2)
  isdublicate=ee.String(ee.Algorithms.If({'condition': cond,
                  'TrueCase': "dublicate",
                  'FalseCase': "no_dublicate"}))
    return image.set({"status_dublicate": isdublicate})


imgcoll_added=imgcollection.map(detect_dublicates)
imgcollfiltered=imgcoll_added.filter(ee.Filter.eq("status_dublicate", "no_dublicate"))
return imgcollfiltered


def noimage():
        #  Create Empty ImageCollection
            emptyimgcol=ee.ImageCollection(ee.Image(5)).filter(ee.Filter.eq('system:index', 0))
  return emptyimgcol


ifNone=ee.Algorithms.If(condition1, normal(), noimage())
return(ifNone)

# merge does not wirk with at the moment! (find out, how an empty imagecollection can be merged)
landsat5datafi=ee.ImageCollection(startdublicatfiltering(landsat5datapre))
landsat7datafi=ee.ImageCollection(startdublicatfiltering(landsat7datapre))
landsat8datafi=ee.ImageCollection(startdublicatfiltering(landsat8datapre))
sentineldatafi=ee.ImageCollection(startdublicatfiltering(sentineldatapre))

# Merge all individual collections per Sensor to one collection
imagecol_all=landsat5datafi.merge(landsat7datafi).merge(landsat8datafi).merge(sentineldatafi)

return imagecol_all
