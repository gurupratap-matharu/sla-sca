/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var sentinel = ee.ImageCollection("COPERNICUS/S2"),
    RGI = ee.FeatureCollection("users/josiaszeller/RGVI_v6");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
// Image-Availability

var glimsid = 'G007871E45994N'
var geometry = RGI.filterMetadata('GLIMSId','equals',glimsid);
var ImageCollection = sentinel.filterBounds(geometry).filterDate('2020-11-27', '2020-12-01')
print(ImageCollection)


// Hello