/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var l8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_TOA"),
  l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_TOA"),
  l5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_TOA"),
  s2 = ee.ImageCollection("COPERNICUS/S2"),
  ING = ee.FeatureCollection(
    "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve"
  );
/***** End of imports. If edited, may not auto-convert in the playground. *****/

exports.imageCollection = function (
  ing_id,
  startyear,
  endyear,
  cloudiness,
  coverage,
  filterDOYstart,
  filterDOYend,
  hsboolean,
  dem
) {
  //////////////////    RENAME BANDS    /////////////////////////////////////////
  var s2_bands = ee.List([
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B8A",
    "B9",
    "B10",
    "B11",
    "B12",
    "QA60",
  ]);
  var s2_band_names = ee.List([
    "cb",
    "blue",
    "green",
    "red",
    "re1",
    "re2",
    "re3",
    "nir",
    "re4",
    "vapor",
    "cirrus",
    "swir1",
    "swir2",
    "BQA",
  ]);

  var l8_bands = ee.List([
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B10",
    "B11",
    "BQA",
  ]);
  var l8_band_names = ee.List([
    "cb",
    "blue",
    "green",
    "red",
    "nir",
    "swir1",
    "swir2",
    "tir1",
    "tir2",
    "BQA",
  ]);

  var l7_bands = ee.List([
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6_VCID_1",
    "B7",
    "BQA",
  ]);
  var l7_band_names = ee.List([
    "blue",
    "green",
    "red",
    "nir",
    "swir1",
    "tir1",
    "swir2",
    "BQA",
  ]);

  var l5_bands = ee.List(["B1", "B2", "B3", "B4", "B5", "B6", "B7", "BQA"]);
  var l5_band_names = ee.List([
    "blue",
    "green",
    "red",
    "nir",
    "swir1",
    "tir1",
    "swir2",
    "BQA",
  ]);
  ///////////////////   PREPARE AND FILTER RGI INVENTORY V6   /////////////////////

  var geometry = ING.filterMetadata("ID_local", "equals", ing_id); //filter the ING inventory with the input of the ID_local --> get the geometry

  //////////////////    PREPARE FILTERING  OF SAT-DATASETS
  var startdate = ee.Date.fromYMD(startyear, 1, 1); //first day of the year
  var enddate = ee.Date.fromYMD(endyear, 12, 31); //last day of the year

  var geometryraw = geometry.geometry();
  var geometrybuffered = geometryraw.buffer(2500, 5);

  ////////////////////////////CLOUD SCORE FOR MASKING///////////////////////////////////////////
  // Compute a cloud score.  This expects the input image to have the common band names
  // Cloud Score is only for the nocloud-ratio and has no influence on the classification. If the cloudratio detected with this score to high,
  // it's possible that the images are filtered within the next step. But: All images are included in the collection
  var cloudScore = function (img) {
    // A helper to apply an expression and linearly rescale the output.
    var rescale = function (img, exp, thresholds) {
      return img
        .expression(exp, { img: img })
        .subtract(thresholds[0])
        .divide(thresholds[1] - thresholds[0]);
    };

    // Compute several indicators of cloudyness and take the minimum of them.
    var score = ee.Image(1.0);
    // Clouds are reasonably bright in the blue band.
    score = score.min(rescale(img, "img.blue", [0.1, 0.3]));

    // Clouds are reasonably bright in all visible bands.
    score = score.min(
      rescale(img, "img.red + img.green + img.blue", [0.2, 0.8])
    );

    // Clouds are reasonably bright in all infrared bands.
    score = score.min(
      rescale(img, "img.nir + img.swir1 + img.swir2", [0.3, 0.8])
    );

    // However, clouds are not snow.
    var ndsi = img.normalizedDifference(["green", "swir1"]);
    var score2 = score.min(rescale(ndsi, "img", [0.7, 0.6]));

    var score3 = ee.Image(1).subtract(score2).select([0], ["cloudscore"]);

    var img1_1 = img.addBands(score3); //add band with cloud score per pixel
    var img1_2 = img1_1.select("cloudscore").gt(0.7); //select cloud pixels to mask original image
    var img1_3 = img.mask(img1_2); //mask image
    return img1_3;
  };
  //function to calculate the area of an image (area of a band without masked pixels)
  var areacalcrast = function (img) {
    var img1 = img.select("green").neq(0); //select not masked pixels
    var fc1 = img1.reduceToVectors({
      //produces a FeatureCollection
      reducer: ee.Reducer.countEvery(),
      geometry: geometry,
      scale: 30,
      maxPixels: 1e10,
      bestEffort: true,
      crs: "EPSG:4326", //use EPSG:4326 -_> explicit that thers no confusion and no check needed for all images
    });
    var area = fc1.reduceColumns({
      reducer: ee.Reducer.sum(),
      selectors: ["count"],
    });
    var area1 = ee.Number(area.get("sum")).divide(1000); //return the number of the area
    return area1; //returns a feature
  };

  //////////////////////////////    FUNCTION: CALCULATE AREA    ///////////////////////////////
  var areacalc = function (image4) {
    //function to calculate the area of an image (area of a band without masked pixels)
    var prep = image4.select("green"); //take one bands (here green) to calculate the area.
    var pixelarea = prep.multiply(ee.Image.pixelArea());

    var areacalc = pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 30,
      bestEffort: true,
      maxPixels: 1e10,
    });
    var numb = areacalc.getNumber("green");
    return numb; //return the number of the area
  };

  ////////////////////////////   CALCULATE GLACIER AREA   //////////////////////////////////////
  var areatotal1 = ee
    .Image(1)
    .clip(geometry)
    .select("constant")
    .rename("green");

  var areatotal = geometry.first().get("Area"); //areacalc(areatotal1);
  var areaIMGglacier1 = areacalcrast(areatotal1); //Var is a number that is included in the prop. for each image

  //////////////////////      Function to add ALBEDO to every image   ////////////////
  //Based on LIANG (2000) - Narrowband to broadband conversions of land surface albedo: I Algorithms
  var addAlbedo = function (img) {
    var albedo = img
      .select("blue")
      .multiply(0.356)
      .add(img.select("red").multiply(0.13))
      .add(img.select("nir").multiply(0.373))
      .add(img.select("swir1").multiply(0.085))
      .add(img.select("swir2").multiply(0.072))
      .subtract(0.0018)
      .rename(["albedo"])
      .divide(1.016);
    var imgexp = img.addBands(albedo);
    return imgexp;
  };

  ////////////////////////////FUNCTION: ADD CLOUD-NoCLOUD-RATIO PER IMAGE////////////////////////
  var cloudfilterlandsat = function (data) {
    var input = data.clip(geometry); //get the right geometry
    var cloudfreefunc = cloudScore(input); //calculate the score
    var countcf = areacalc(cloudfreefunc); //calculate the area for the image when the detected clouds are masked
    var countnorm = areacalc(input); //calculate the area for the image (with clouds)
    var nocloud = countcf.divide(countnorm).multiply(100); //calculate ratio
    return nocloud; //return number
  };

  //function to add all the 'Quality-Information' to the image
  var addquainfo = function (data) {
    //data is one image from the imagecollection (unfiltered)
    var clip = data.clip(geometrybuffered); //clip image
    var quality = cloudfilterlandsat(clip); //calculate the cloudscore with the function above, save number with the ratio
    var image = clip.set("noclouds", quality); //set cloudratio to the image properties
    var Sens = clip.get("SPACECRAFT_ID");
    var areaRGIoutline = areacalcrast(clip); //calculate the area for the image to detect, if there are parts without information
    var areaIMGglacier = areaIMGglacier1; //areacalcrast(areatotal1)
    var arearatio = areaRGIoutline.divide(areaIMGglacier).multiply(100); //calculate ration of coverage from one image over the glacier
    var image1 = image
      .set("ING_ID", ing_id) //set ING to image
      .set("areaRGIoutline", areaRGIoutline)
      .set("areaIMGglacier", areaIMGglacier)
      .set("arearatio", arearatio)
      .set("SENSOR", Sens)
      .set("hsboolean", hsboolean)
      .set("deminfo", dem);

    return image1;
  };

  var cloudfiltersentinel = function (data) {
    var input = data.clip(geometry);
    var cloudfreefunc = cloudScore(input);
    var countcf = areacalc(cloudfreefunc);
    var countnorm = areacalc(input);
    var nocloud = countcf.divide(countnorm).multiply(100);
    return nocloud;
  };

  var saddquainfo = function (data) {
    var clip = data.clip(geometrybuffered);
    var quality = cloudfiltersentinel(clip);
    var areaINGoutline = areacalcrast(clip);
    var areaIMGglacier = areaIMGglacier1;
    var arearatio = areaINGoutline.divide(areaIMGglacier).multiply(100);
    var rightangle = ee.Number(90);
    var sunAZ = ee.Number(clip.get("MEAN_SOLAR_AZIMUTH_ANGLE"));
    var sunELZenith = ee.Number(clip.get("MEAN_SOLAR_ZENITH_ANGLE"));
    var Sens = clip.get("SPACECRAFT_NAME");
    var sunEL = rightangle.subtract(ee.Number(sunELZenith));
    var image = clip.set("noclouds", quality);
    var image1 = image
      .set("ING_ID", ing_id)
      .set("areaINGoutline", areaINGoutline)
      .set("areaIMGglacier", areaIMGglacier)
      .set("arearatio", arearatio)
      .set("SUN_AZIMUTH", sunAZ)
      .set("SUN_ELEVATION", sunEL)
      .set("SENSOR", Sens)
      .set("hsboolean", hsboolean)
      .set("deminfo", dem);
    return image1;
  };

  ////////////////    RUN ALL FILTERS FOR EVERY SENSOR    ///////////////////////////////////
  var cloudfilter1 = ee.Filter.gt("noclouds", cloudiness); //get the threshold (in %) from the input and create the cloudfilter
  var coveragefilter = ee.Filter.gt("arearatio", coverage); //get threshold in % from input

  var landsat5data = l5
    .filterBounds(geometry) //filter on the base of the choosen geometry
    .filterDate(startdate, enddate)
    .filter(
      ee.Filter.calendarRange(filterDOYstart, filterDOYend, "day_of_year")
    )
    .map(function (img) {
      //create image collection
      var img1 = img.select(l5_bands).rename(l5_band_names); //rename all images with the correct names
      var img2 = img1; //.reproject('EPSG:32632',null, 30);
      return img2;
    })
    .map(addquainfo) //map all quality information and set it as properties for each image
    .filter(cloudfilter1) //filter the images with the cloud-threshold
    .filter(coveragefilter)
    .map(addAlbedo);

  var landsat7data = l7
    .filterBounds(geometry)
    .filterDate(startdate, enddate)
    .filter(
      ee.Filter.calendarRange(filterDOYstart, filterDOYend, "day_of_year")
    )
    .map(function (img) {
      var img1 = img.select(l7_bands).rename(l7_band_names);
      var img2 = img1; //.reproject('EPSG:32632',null, 30);
      return img2;
    })
    .map(addquainfo)
    .filter(cloudfilter1)
    .filter(coveragefilter)
    .map(addAlbedo);

  var landsat8data = l8
    .filterBounds(geometry)
    .filterDate(startdate, enddate)
    .filter(
      ee.Filter.calendarRange(filterDOYstart, filterDOYend, "day_of_year")
    )
    .map(function (img) {
      var img1 = img.select(l8_bands).rename(l8_band_names);
      var img2 = img1; //.reproject('EPSG:32632',null, 30);
      return img2;
    })
    .map(addquainfo)
    .filter(cloudfilter1)
    .filter(coveragefilter)
    .map(addAlbedo);

  var sentineldata = s2
    .filterBounds(geometry)
    .filterDate(startdate, enddate)
    .filter(
      ee.Filter.calendarRange(filterDOYstart, filterDOYend, "day_of_year")
    )
    .map(function (img) {
      var img1 = img.select(s2_bands).rename(s2_band_names);
      var img2 = img1; //.reproject('EPSG:32632',null, 30);
      return img2;
    })

    .map(saddquainfo)
    .filter(cloudfilter1)
    .filter(coveragefilter)
    .map(function (img) {
      //scale the pixelvalues to a ratio between 0 and 1
      var img1 = img
        .select([
          "cb",
          "blue",
          "green",
          "red",
          "re1",
          "re2",
          "re3",
          "nir",
          "re4",
          "vapor",
          "cirrus",
          "swir1",
          "swir2",
        ])
        .divide(10000);
      return img.addBands(img1, null, true);
    })
    .map(addAlbedo);

  //////////////////////   FILTER DUBLICATES WITH SAME DATE    //////////////////////////////////////////////////
  //sort all collection (for each sensor), that the dublicate filter works
  var landsat5datapre = landsat5data.sort("system:time_start");
  var landsat7datapre = landsat7data.sort("system:time_start");
  var landsat8datapre = landsat8data.sort("system:time_start");
  var sentineldatapre = sentineldata.sort("system:time_start");
  // Function to detect dublicates
  var startdublicatfiltering = function (imgcollection) {
    var condition1 = imgcollection.size();

    var normal = function () {
      var list = imgcollection.toList(imgcollection.size());
      var image = ee.Image(list.get(0));
      var dummyimage = ee.Image(1).set("system:time_start", 0);
      //Add in the end of the list a dummy image
      list = list.add(dummyimage);

      var detect_dublicates = function (image) {
        var isdublicate = ee.String("");
        var number = list.indexOf(image);
        var image1 = ee.Image(list.get(number.add(1)));
        //Compare the image(0) in the ImageCollection with the image(1) in the List
        var date1 = image.date().format("Y-M-d");
        var date2 = image1.date().format("Y-M-d");
        var cond = ee.Algorithms.IsEqual(date1, date2);
        isdublicate = ee.String(
          ee.Algorithms.If({
            condition: cond,
            trueCase: "dublicate",
            falseCase: "no_dublicate",
          })
        );
        return image.set({ status_dublicate: isdublicate });
      };

      var imgcoll_added = imgcollection.map(detect_dublicates);
      var imgcollfiltered = imgcoll_added.filter(
        ee.Filter.eq("status_dublicate", "no_dublicate")
      );
      return imgcollfiltered;
    };
    var noimage = function () {
      //  Create Empty ImageCollection
      var emptyimgcol = ee
        .ImageCollection(ee.Image(5))
        .filter(ee.Filter.eq("system:index", 0));
      return emptyimgcol;
    };

    var ifnull = ee.Algorithms.If(condition1, normal(), noimage());
    return ifnull;
  };
  //merge does not work with at the moment! (find out, how an empty imagecollection can be merged)
  var landsat5datafi = ee.ImageCollection(
    startdublicatfiltering(landsat5datapre)
  );
  var landsat7datafi = ee.ImageCollection(
    startdublicatfiltering(landsat7datapre)
  );
  var landsat8datafi = ee.ImageCollection(
    startdublicatfiltering(landsat8datapre)
  );
  var sentineldatafi = ee.ImageCollection(
    startdublicatfiltering(sentineldatapre)
  );

  // Merge all individual collections per Sensor to one collection
  var imagecol_all = landsat5datafi
    .merge(landsat7datafi)
    .merge(landsat8datafi)
    .merge(sentineldatafi);

  return imagecol_all;
};

// hello
