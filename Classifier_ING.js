/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var ALOS = ee.Image("JAXA/ALOS/AW3D30/V2_2"),
    srtm = ee.Image("USGS/SRTMGL1_003"),
    ING = ee.FeatureCollection("users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
 
exports.decisiontree = function(img) {  
 
//// get input values from the image to be classified
var hsboolean =   img.get('hsboolean'); //check boolean (added to every image as a property)
var ing_id =     img.get('ING_ID') ; //get GLIMS-ID
var deminfo =     img.get('deminfo');
var geometry =    ING.filterMetadata('ID_local','equals',ing_id); //Get Shapefile


//------------DEM-Selection
var dem1 = srtm.select('elevation').rename('AVE_DSM').clip(geometry); //rename srtm dem to match the alos dem, clip to geometry
var dem2 = ALOS.select('AVE_DSM').clip(geometry);
var demselector = ee.Algorithms.IsEqual(ee.String(deminfo), ee.String('ALOS')); //boolean to check, if dem selection == ALOS
var demselection = ee.Image(ee.Algorithms.If(demselector, dem2, dem1)); //select the DEm that should be used
var demglacier = demselection; //assignment
//------------

//cast geometry
var geometryraw = geometry.geometry();
//buffer geometry withing 2500 m to clip the DSM
var geometrybuffered = geometryraw.buffer(2500,5);
var elevation=demglacier.clip(geometrybuffered);
//get the CRS-Information from the original image (Landsat or Sentinel-2)
var crstransform = img.select('blue').projection();



var ndsi1 = img.normalizedDifference(['green', 'swir1']);
var ndwi1 = img.normalizedDifference(['green', 'nir']);
 
var ndsi = ndsi1.select('nd').rename('ndsi');
var ndwi = ndwi1.select('nd').rename('ndwi');

var img2 = img.addBands(ndsi);
var img3 = img2.addBands(ndwi);

///////   Hill Shadow   ///////////////////////////////////////////////
// check, if hillshadow-map should added or not
var Yaddhillshadow = function() { //function to add a layer with the hillshadow
var rightangle = ee.Number(90); //add 90 to calculate difference
var zenith = rightangle.subtract(ee.Number(img.get('SUN_ELEVATION'))); //just a calculation to get the zenit-value in degree
var shadowMap=ee.Terrain.hillShadow(elevation, img.get('SUN_AZIMUTH'), zenith,200,true)
                                        .focal_min(4).reproject(crstransform);
var img4=img3.addBands(shadowMap); //add band (shadow band) to the original image
return img4;
};

var Naddhillshadow = function(){ //function to add a dummy-layer with <-1> value for every pixel
var shadowMap = ee.Image(-1).clip(geometry).rename(['shadow']);//.reproject(crstransform); //create dummy image with minus 1 values
var img4=img3.addBands(shadowMap); //add band (shadow band) to the original image
return img4;
};
// function to select the right function (add hill shadow band or not, depending on the properties of the image
var addhillshadow = ee.Algorithms.If(hsboolean,Yaddhillshadow(),Naddhillshadow());
var img4 = ee.Image(addhillshadow);
var systemtime =  ee.Number(img.get('system:time_start'));
var hillshadowmask =  img4
                     .select('shadow')
                     .eq(0)
                     .clip(geometry);

///////   Start Classification   ///////////////////////////////////////////////

// Blueprint
var blueprint = ee.Image(1).select('constant').rename('classification');
var start = img4.addBands(blueprint);

//get all pixels, that are marked as cloudshadows from previous step
var cloudshadowmask = img.select('cloudShadow').eq(1);
//combine cloudshadows and hillshadow into one mask
var shadowmaskhillandcloud = cloudshadowmask.add(hillshadowmask).gte(1);
var hillshadow1 = start.mask(shadowmaskhillandcloud); //
var rescale2 = function(img, exp, thresholds) {
               return img.expression(exp, {img: img})
                      .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
                      };
///////   Create Shadowmask (Score)   ///////////////////////////////////////////////
  var score2 = ee.Image(1.0);
// Snow is reasonably bright in all visible bands.
  score2 = score2.min(rescale2(hillshadow1, 'img.red + img.green + img.blue', [0.5, 2]));

var shadowonsnowmask = score2.select('constant').gt(0.45);
var undefshadowmask = score2.select('constant').lt(0.45);

var shadowonsnow1 = hillshadow1.mask(shadowonsnowmask).select('classification').multiply(6);
var shadowonsnow = shadowonsnow1.unmask(-1).clip(geometry);
var undefshadow1 = hillshadow1.mask(undefshadowmask).select('classification').multiply(8);
var undefshadow = undefshadow1.unmask(-1).clip(geometry);

//Prepare image without shaded area (shadow is masked)
var img5 = start.mask((shadowmaskhillandcloud).not());

var watermask = img5.select('ndwi').gt(0.3); //create watermask with NDWI-Values higher than 0.3
var notwater = img5.select('ndwi').lt(0.3); //creat mask to carry on the classification with everything exept the water
var waterclass = start.mask(watermask);

var shadowmask = waterclass.select('blue').gt(0.1);
var shadow1 = start.mask(shadowmask).select('classification').multiply(8);
var shadow = shadow1.unmask(-1); //mask everything else with minus 1 (not a class)

var watermask1 = waterclass.select('blue').lt(0.1);
var water1 = start.mask(watermask1).select('classification').multiply(2);
var water = water1.unmask(-1); //mask everything else with minus 1 (not a class)

var img55 = img4.mask(notwater);

//Differentiate between snow/ice and clouds/debris

var snowice = img55.select('ndsi').gt(0.4);
var clouddebris = img55.select('ndsi').lt(0.4);

var img6 = img55.mask(snowice);
///////   OTSU THRESHOLD FOR NIR-BAND   ///////////////////////////////////////////////

// Compute the histogram of the NIR band.

var otsuimg = img6.clip(geometry).select('nir').multiply(10000);
var histogram = otsuimg.reduceRegion({
  reducer: ee.Reducer.histogram(50,2)
      .combine('mean', null, true)
      .combine('variance', null, true), 
  geometry: geometry, 
  scale: 30,
  bestEffort: true
});

// Return the DN that maximizes interclass variance in NIR
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  
  // Compute between sum of squares, where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
// Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};
var nirthreshold = otsu(histogram.get('nir_histogram')).divide(10000); 

var otsu1 = nirthreshold   ;                                  // get threshold to check if it is to0 low or too high
var checklow = otsu1.gt(0.41);                                // set lower boundary
var checkhigh = otsu1.lt(0.54);                               // set higher boundary
var checkotsu = checklow.add(checkhigh);                      // add the two boolean. Between the two boundaries if ==2, otherwise >2
var useotsu = checkotsu.gt(1);                                // check if otsu between boundaries

var fixed = 0.47                         ;                    // fixed threshold if otsu out of range
var threshold = ee.Algorithms.If(useotsu,otsu1,fixed);        // write threshold: fixed or use the original otsu-number
var usethis = ee.Number(threshold);

///////   Mask Snow and Ice   ///////////////////////////////////////////////
var snowmask = img6.select('nir').gt(usethis); 
var snowclass = start.mask(snowmask).select('classification').multiply(1);
var snow = snowclass.unmask(-1);

var icemask = img6.select('nir').lt(usethis);

// reclassifie dark snow to ice
var icesnowdif = start.mask(icemask);
var rescale = function(img, exp, thresholds) {
               return img.expression(exp, {img: img})
                      .subtract(thresholds[0]).divide(thresholds[1] - thresholds[0]);
                      };
  var score = ee.Image(1.0);
// Snow is reasonably bright in all visible bands.
  score = score.min(rescale(icesnowdif, 'img.red + img.green + img.blue', [0.5, 2]));

var snow2 = score.select('constant').gt(0.45);
var icemask2 = score.select('constant').lt(0.45);
var iceclass = start.mask(icemask2).select('classification').multiply(0);
var ice = iceclass.unmask(1);




///////   Mask Cloud and Debriscover   ///////////////////////////////////////////////

var img7 = img5.mask(clouddebris);
 
var cloudmask = img7.select('red').gt(0.6);
var cloudclass = start.mask(cloudmask).select('classification').multiply(4);
var cloud = cloudclass.unmask(-1);

var debrismask = img7.select('red').lt(0.3);
var debrisclass = start.mask(debrismask).select('classification').multiply(3);
var debris = debrisclass.unmask(-1);

///////   Reprojection for all Masks   ///////////////////////////////////////////////
var water11 = water.reproject(crstransform);
var cloud1 = cloud.reproject(crstransform);
var debris1 = debris.reproject(crstransform);
var shadow111 = shadow.reproject(crstransform);
var shadowonsnow111 = shadowonsnow.reproject(crstransform);
var snow1 = snow.reproject(crstransform);
var undefshadow111 = undefshadow.reproject(crstransform);
var ice1 = ice.reproject(crstransform);

///////   Create one single image with all masks   ///////////////////////////////////////////////
var classify = water11.max(cloud1).max(debris1)
                    .max(shadow111).max(shadowonsnow111)
                    .max(snow1).max(undefshadow111)
                    .max(ice1)//.reproject(crstransform);
                    .clip(geometry);

var vectorizedmap = classify.reduceToVectors({
  scale:            30,
  geometryType:     'polygon',
  geometry:         geometry,
  reducer:          ee.Reducer.countEvery(),
  bestEffort:       true
});


/////////////////   START RF-CLASSSIFICATION FOR ICE/SNOW-AREA    ////////////////////

var s2boolean = ee.Algorithms.IsEqual(img.get('SENSOR'),'Sentinel-2A');
var l5boolean = ee.Algorithms.IsEqual(img.get('SENSOR'),'LANDSAT_5');
var l7boolean = ee.Algorithms.IsEqual(img.get('SENSOR'),'LANDSAT_7');
var l8boolean = ee.Algorithms.IsEqual(img.get('SENSOR'),'LANDSAT_8');

var s2bands = ee.List(['green', 'blue','red','nir','swir1', 'swir2', 'albedo']);
var l7bands = ee.List(['green', 'blue','red','nir','swir1', 'swir2', 'albedo']);
var l8bands = ee.List(['green', 'blue','red','nir','swir1', 'swir2', 'albedo']);
var l5bands = ee.List(['green', 'blue','red','nir','swir1', 'swir2', 'albedo']);
var bandsbackup = ee.List(['green', 'blue','red','nir','swir1', 'swir2', 'albedo']);

var l8check = ee.Algorithms.If(l8boolean, l8bands, bandsbackup);
var l7check = ee.Algorithms.If(l7boolean, l7bands, l8check);
var l5check = ee.Algorithms.If(l5boolean, l5bands, l7check);
var bnds = ee.Algorithms.If(s2boolean, s2bands, l5check);


//select only ice polygons
var icevec = vectorizedmap.filter(ee.Filter.eq('label', 0));
//select only snow polygons
var snowvec = vectorizedmap.filter(ee.Filter.eq('label', 1));
var watervec = vectorizedmap.filter(ee.Filter.eq('label', 2));
var debrisvec = vectorizedmap.filter(ee.Filter.eq('label', 3));
var cloudsvec = vectorizedmap.filter(ee.Filter.eq('label', 4));
var shadsnowvec = vectorizedmap.filter(ee.Filter.eq('label', 6));
var shadowvec = vectorizedmap.filter(ee.Filter.eq('label', 8));

//set <<points>> in the vectorized areas for all available classes.
var points = 100
var randomPointsice = ee.FeatureCollection
          .randomPoints(icevec, points)
          .map(function(feat){
                                return feat.set('cover', 0);
                                  });
var randomPointssnow = ee.FeatureCollection
          .randomPoints(snowvec, points)
          .map(function(feat){
                                return feat.set('cover', 1);
                                  });
var randomPointsdebris = ee.FeatureCollection
          .randomPoints(debrisvec, points)
          .map(function(feat){
                                return feat.set('cover', 3);
                                  });
var randomPointswater = ee.FeatureCollection
          .randomPoints(watervec, points)
          .map(function(feat){
                                return feat.set('cover', 2);
                                  });
var randomPointsclouds = ee.FeatureCollection
          .randomPoints(cloudsvec, points)
          .map(function(feat){
                                return feat.set('cover', 4);
                                  });
var randomPointshadsnow = ee.FeatureCollection
          .randomPoints(shadsnowvec, points)
          .map(function(feat){
                                return feat.set('cover', 6);
                                  });
var randomPointshadow = ee.FeatureCollection
          .randomPoints(shadowvec, points)
          .map(function(feat){
                                return feat.set('cover', 8);
                                  });
                                
//merge all random set points in one Feature Collection
var fc1 = randomPointsice.merge(randomPointssnow).merge(randomPointsdebris).merge(randomPointswater).merge(randomPointsclouds)
.merge(randomPointshadsnow).merge(randomPointshadow);
var classProperty = 'cover';
var bands = bnds;


// Prepare and Train the classifier.
var training1_data = img.select(bands).sampleRegions({
  collection: fc1,
  properties: null,
  scale: 30 });


var classifier = ee.Classifier.smileRandomForest(10).train({
  features: training1_data,
  classProperty: classProperty,
  inputProperties: bands,

});
//start classification with the trained random forest-classifier
var classifiedRF = img.clip(geometry).classify(classifier);
//.reproject(crstransform);
 
var metadatafc = function(img){
        var properties = img.propertyNames();
        var addMeta = function(prop){
          var feat = ee.Feature(null);
          var name = prop;
          var value = img.get(prop);
          var meta = feat.set(name, value);
          return(meta);
        };
        var META = properties.map(addMeta);
        return(META);
};
var fcmetadata = metadatafc(img);

var classified = classifiedRF.set('system_time_start',systemtime)
                          .set('ING_ID', ing_id).set('otsu',nirthreshold).set('deminfo', deminfo)
                          .clip(geometry).set('metadata', fcmetadata)
                          .mask((img.select('blue').gte(0)));
                         // .reproject(crstransform);
return classified;
};

// Hello