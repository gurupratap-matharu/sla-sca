"""
Preparation Module for SLA - SCA

This script composes of various helper functions but its main role is to prepare and provide an image collection. The script exposes a high level
function which does all the heavy lifting.

The composite image is generated with this script and used in main code for further processing.
"""

import ee
from ee_plugin import Map

ING = ee.FeatureCollection(
    "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve")
L5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA")
L7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_TOA")
L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA")
S2 = ee.ImageCollection("COPERNICUS/S2")


def imageCollection(ing_id, startyear, endyear, cloudiness, coverage, filterDOYstart, filterDOYend, hsboolean, dem):

    l5_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'BQA'])
    l5_band_names = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'tir1', 'swir2', 'BQA'])

    l7_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6_VCID_1', 'B7', 'BQA'])
    l7_band_names = ee.List(['blue', 'green', 'red', 'nir', 'swir1', 'tir1', 'swir2', 'BQA'])

    l8_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B10', 'B11', 'BQA'])
    l8_band_names = ee.List(['cb', 'blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'tir1', 'tir2', 'BQA'])

    s2_bands = ee.List(['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12', 'QA60'])
    s2_band_names = ee.List(['cb', 'blue', 'green', 'red', 're1', 're2', 're3',
                            'nir', 're4', 'vapor', 'cirrus', 'swir1', 'swir2', 'BQA'])

    geometry = ING.filterMetadata('ID_local', 'equals', ing_id)

    start_date = ee.Date.fromYMD(startyear, 1, 1)
    end_date = ee.Date.fromYMD(endyear, 12, 31)

    geometry_raw = geometry.geometry()
    geometry_buffered = geometry_raw.buffer(2500, 5)


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


# function to calculate the area of an image (area of a band without masked pixels)
def areacalcrast(img):
                img1 = img.select('green').neq(0);  # select not masked pixels
                fc1 = img1.reduceToVectors({  # produces a FeatureCollection
                                      'reducer': ee.Reducer.countEvery(),
                                      'geometry': geometry,
                                      'scale': 30,
                                      'maxPixels': 1e10,
                                      'bestEffort': True,
                                      'crs': 'EPSG:4326'  # use EPSG:4326 -_> explicit that thers no confusion and no check needed for all images
                                  })
               area=fc1.reduceColumns({
                                      'reducer': ee.Reducer.sum(),
                                      'selectors': ['count'],
                                  })
              area1=ee.Number(area.get('sum')).divide(1000);  # return the number of the area
            return area1;  # returns a feature



# FUNCTION: CALCULATE AREA    ###############/
def areacalc(image4)  # to calculate the area of an image (area of a band without masked pixels):
                prep=image4.select('green');  # take one bands (here green) to calculate the area.
                pixelarea=prep.multiply(ee.Image.pixelArea())

                  areacalc=pixelarea.reduceRegion({
                       'reducer': ee.Reducer.sum(),
                       'geometry': geometry,
                       'scale': 30,
                       'bestEffort': True,
                       'maxPixels': 1e10
                  })
              numb=areacalc.getNumber('green')
              return numb;  # return the number of the area


##############   CALCULATE GLACIER AREA   ###################
areatotal1=ee.Image(1).clip(geometry).select('constant').rename('green')

areatotal=geometry.first().get('Area');  # areacalc(areatotal1)
areaIMGglacier1=areacalcrast(areatotal1);  # Var is a number that is included in the prop. for each image

###########      Function to add ALBEDO to every image   ########
# Based on LIANG (2000) - Narrowband to broadband conversions of land surface albedo: I Algorithms
def addAlbedo(img):
  albedo=(((img.select('blue').multiply(0.356)) \
                .add(img.select('red').multiply(0.130)) \
                .add(img.select('nir').multiply(0.373)) \
                .add(img.select('swir1').multiply(0.085)) \
                .add(img.select('swir2').multiply(0.072)) \
                .subtract(0.0018).rename(['albedo'])).divide(1.016))
  imgexp=img.addBands(albedo)
  return imgexp


##############FUNCTION: ADD CLOUD-NoCLOUD-RATIO PER IMAGE############
def cloudfilterlandsat(data):
                    input=data.clip(geometry);  # get the right geometry
                    cloudfreefunc=cloudScore(input);  # calculate the score
                    # calculate the area for the image when the detected clouds are masked
                    countcf=areacalc(cloudfreefunc);
                    countnorm=areacalc(input);  # calculate the area for the image (with clouds)
                    nocloud=countcf.divide(countnorm).multiply(100);  # calculate ratio
                return nocloud;  # return number


# function to add all the 'Quality-Information' to the image
def addquainfo(data)  # data is one image from the imagecollection (unfiltered):
                    clip=data.clip(geometrybuffered);  # clip image
                    # calculate the cloudscore with the function above, save number with the ratio
                    quality=cloudfilterlandsat(clip);
                    image=clip.set('noclouds', quality);  # set cloudratio to the image properties
                    Sens=clip.get('SPACECRAFT_ID')
                    # calculate the area for the image to detect, if there are parts without information
                    areaRGIoutline=areacalcrast(clip);
                    areaIMGglacier=areaIMGglacier1;  # areacalcrast(areatotal1)
                    # calculate ration of coverage from one image over the glacier
                    arearatio=areaRGIoutline.divide(areaIMGglacier).multiply(100);
                    image1=image  .set('GLIMSID', glimsid)  # set glimsid to image \
                                        .set('areaRGIoutline', areaRGIoutline) \
                                        .set('areaIMGglacier', areaIMGglacier) \
                                        .set('arearatio', arearatio) \
                                        .set('SENSOR', Sens) \
                                        .set('hsboolean', hsboolean) \
                                        .set('deminfo', dem)

              return image1


def cloudfiltersentinel(data):
                    input=data.clip(geometry)
                    cloudfreefunc=cloudScore(input)
                    countcf=areacalc(cloudfreefunc)
                    countnorm=areacalc(input)
                    nocloud=countcf.divide(countnorm).multiply(100)
              return nocloud


def saddquainfo(data):
                    clip=data.clip(geometrybuffered)
                    quality=cloudfiltersentinel(clip)
                    areaRGIoutline=areacalcrast(clip)
                    areaIMGglacier=areaIMGglacier1
                    arearatio=areaRGIoutline.divide(areaIMGglacier).multiply(100)
                    rightangle=ee.Number(90)
                    sunAZ=ee.Number(clip.get('MEAN_SOLAR_AZIMUTH_ANGLE'))
                    sunELZenith=ee.Number(clip.get('MEAN_SOLAR_ZENITH_ANGLE'))
                    Sens=(clip.get('SPACECRAFT_NAME'))
                    sunEL=rightangle.subtract(ee.Number(sunELZenith))
                    image=clip.set('noclouds', quality)
                    image1=image.set('GLIMSID', glimsid) \
                        .set('areaRGIoutline', areaRGIoutline) \
                        .set('areaIMGglacier', areaIMGglacier) \
                        .set('arearatio', arearatio) \
                        .set('SUN_AZIMUTH', sunAZ) \
                        .set('SUN_ELEVATION', sunEL) \
                        .set('SENSOR', Sens) \
                        .set('hsboolean', hsboolean) \
                        .set('deminfo', dem)
              return image1



# RUN ALL FILTERS FOR EVERY SENSOR    #################/
cloudfilter1=ee.Filter.gt('noclouds', cloudiness);  # get the threshold (in %) from the input and create the cloudfilter
coveragefilter=ee.Filter.gt('arearatio', coverage);  # get threshold in % from input

landsat5data=l5.filterBounds(geometry)  # filter on the base of the choosen geometry \
                .filterDate(startdate, enddate) \
                .filter(ee.Filter.calendarRange(filterDOYstart, filterDOYend, 'day_of_year'))

def func_owr(img)  # create image collection:
                  img1=img.select(l5_bands).rename(l5_band_names);  # rename all images with the correct names
                  img2=img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_owr) \
                .map(addquainfo) \
                .filter(cloudfilter1) \
                .filter(coveragefilter) \
                .map(addAlbedo)


landsat7data=l7.filterBounds(geometry) \
                  .filterDate(startdate, enddate).filter(ee.Filter.calendarRange(filterDOYstart, filterDOYend, 'day_of_year'))

def func_xbh(img):
                  img1=img.select(l7_bands).rename(l7_band_names)
                  img2=img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                  .map(func_xbh) \
                .map(addquainfo) \
                .filter(cloudfilter1) \
                .filter(coveragefilter) \
                .map(addAlbedo)


landsat8data=l8.filterBounds(geometry) \
                .filterDate(startdate, enddate).filter(ee.Filter.calendarRange(filterDOYstart, filterDOYend, 'day_of_year'))

def func_tfu(img):
                  img1=img.select(l8_bands).rename(l8_band_names)
                  img2=img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_tfu) \
                .map(addquainfo) \
                .filter(cloudfilter1) \
                .filter(coveragefilter) \
                .map(addAlbedo)


sentineldata=s2.filterBounds(geometry) \
                .filterDate(startdate, enddate).filter(ee.Filter.calendarRange(filterDOYstart, filterDOYend, 'day_of_year'))

def func_ngk(img):
                  img1=img.select(s2_bands).rename(s2_band_names)
                  img2=img1  # .reproject('EPSG:32632',None, 30)
                  return img2 \
                .map(func_ngk) \
                .map(saddquainfo) \
                .filter(cloudfilter1) \
                .filter(coveragefilter)

def func_jty(img):
                  # scale the pixelvalues to a ratio between 0 and 1
                  img1=img \
                    .select((['cb', 'blue', 'green', 'red', 're1', 're2', 're3', 'nir', 're4', 'vapor', 'cirrus', 'swir1', 'swir2'])) \
                    .divide(10000)
                  return img.addBands(img1, None, True) \
                .map(func_jty
).map(addAlbedo)





).map(addAlbedo)

###########   FILTER DUBLICATES WITH SAME DATE    #########################
# sort all collection (for each sensor), that the dublicate filter works
landsat5datapre = landsat5data.sort('system:time_start')
landsat7datapre = landsat7data.sort('system:time_start')
landsat8datapre = landsat8data.sort('system:time_start')
sentineldatapre = sentineldata.sort('system:time_start')
# Function to detect dublicates
def startdublicatfiltering(imgcollection):


condition1 = imgcollection.size()


def normal():


list = imgcollection.toList(imgcollection.size())
image = ee.Image(list.get(0))
dummyimage = ee.Image(1).set('system:time_start', 0)
# Add in the end of the list a dummy image
list = list.add(dummyimage)


def detect_dublicates(image):
  isdublicate = ee.String("")
  number = list.indexOf(image)
  image1 = ee.Image(list.get(number.add(1)))
  # Compare the image(0) in the ImageCollection with the image(1) in the List
  date1 = image.date().format("Y-M-d")
  date2 = image1.date().format("Y-M-d")
  cond = ee.Algorithms.IsEqual(date1, date2)
  isdublicate = ee.String(ee.Algorithms.If({'condition': cond,
                  'TrueCase': "dublicate",
                  'FalseCase': "no_dublicate"}))
    return image.set({"status_dublicate": isdublicate})


imgcoll_added = imgcollection.map(detect_dublicates)
imgcollfiltered = imgcoll_added.filter(ee.Filter.eq("status_dublicate", "no_dublicate"))
return imgcollfiltered


def noimage():
        #  Create Empty ImageCollection
            emptyimgcol = ee.ImageCollection(ee.Image(5)).filter(ee.Filter.eq('system:index', 0))
  return emptyimgcol


ifNone = ee.Algorithms.If(condition1, normal(), noimage())
return(ifNone)

# merge does not wirk with at the moment! (find out, how an empty imagecollection can be merged)
landsat5datafi = ee.ImageCollection(startdublicatfiltering(landsat5datapre))
landsat7datafi = ee.ImageCollection(startdublicatfiltering(landsat7datapre))
landsat8datafi = ee.ImageCollection(startdublicatfiltering(landsat8datapre))
sentineldatafi = ee.ImageCollection(startdublicatfiltering(sentineldatapre))

# Merge all individual collections per Sensor to one collection
imagecol_all = landsat5datafi.merge(landsat7datafi).merge(landsat8datafi).merge(sentineldatafi)

return imagecol_all



# hello
