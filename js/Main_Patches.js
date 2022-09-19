/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var RGI = ee.FeatureCollection("users/josiaszeller/RGVI_v6"),
  alos = ee.Image("JAXA/ALOS/AW3D30/V2_2"),
  srtm = ee.Image("USGS/SRTMGL1_003");
/***** End of imports. If edited, may not auto-convert in the playground. *****/

exports.sla_extract_patch = function (img) {
  //import image to extract SLA

  ///////   Creating variables   ///////////////////////////////////////////////
  var glimsid = img.get("GLIMSID");
  var otsu = img.get("otsu");
  var deminfo = img.get("deminfo");
  var geometry = RGI.filterMetadata("GLIMSId", "equals", glimsid);

  var classified = img;

  //------------DEM-Selection
  var dem1 = srtm.select("elevation").rename("AVE_DSM").clip(geometry); //rename srtm dem to match the alos dem, clip to geometry
  var dem2 = alos.select("AVE_DSM").clip(geometry);
  var demselector = ee.Algorithms.IsEqual(
    ee.String(deminfo),
    ee.String("ALOS")
  ); //boolean to check, if dem selection == ALOS
  var demselection = ee.Image(ee.Algorithms.If(demselector, dem2, dem1)); //select the DEm that should be used
  var demglacier = demselection; //assignment
  //------------

  ///////   CALCULATE AREA OF IMAGE AND RETURN AREA   ///////////////////////////////////////////////

  var areacalc = function (image4) {
    //function to calculate the area of an image
    //(area of a band without masked pixels)
    var prep = image4.select("AVE_DSM"); //take one bands (here Elevation) to calculate the area.
    var pixelarea = prep.multiply(ee.Image.pixelArea());
    var areacalc1 = pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 20,
      maxPixels: 1e8,
      bestEffort: true,
    });
    var numb = areacalc1.getNumber("AVE_DSM");
    return numb; //return the number of the area
  };
  //////////////    START CODE    ///////////////////////////////////////////////
  //combine classes from classified image
  var dcmask = classified.select("classification").neq(-1); //everything thats not classified
  var mask1 = classified.select("classification").neq(7); //everything thats not shadow on water
  var mask2 = classified.select("classification").neq(4); //everything thats not clouds
  var mask3 = classified.select("classification").neq(3); //everything thats not debris cover
  var mask4 = classified.select("classification").neq(2); //everything thats not water
  var mask5 = classified.select("classification").neq(8);
  var classmask1 = classified.mask(dcmask).clip(geometry);
  var classmask2 = classmask1.updateMask(mask1).clip(geometry);
  var classmask3 = classmask2.updateMask(mask2).clip(geometry);
  var classmask4 = classmask3.updateMask(mask3).clip(geometry);
  var classmask5 = classmask4.updateMask(mask4).clip(geometry);
  var classmask6 = classmask5.updateMask(mask5).clip(geometry);

  ///////   ELEVATION ANALYSIS   ///////////////////////////////////////////////

  var elevation = demglacier.select("AVE_DSM").clip(geometry);

  //combine snow and shadow on snow
  var singleclasssnow = classified.select("classification").eq(1);
  var singleclasssos = classified.select("classification").eq(6); //shadow on snow => classification value ==6
  var snowmask1 = singleclasssnow.max(singleclasssos);
  //combine ice and shadow on ice
  var singleclassice = classified.select("classification").eq(0);
  var singleclasssoi = classified.select("classification").eq(5);
  var icemask1 = singleclassice.max(singleclasssoi);
  //mask ice
  var elevclassice = elevation.mask(icemask1.clip(geometry));
  //mask snow
  var elevclassnow = elevation.mask(snowmask1.clip(geometry));

  // create rasterimage to vectorize
  var forvec = classified
    .mask(snowmask1)
    .select("classification")
    .multiply(0)
    .unmask(-10)
    .clip(geometry);
  var forvec1 = classified
    .mask(icemask1)
    .select("classification")
    .add(10)
    .unmask(-10)
    .clip(geometry);
  var forvecfinal = forvec.max(forvec1);

  ///////   Calculate and Store Areas for Ratio-Calculations   ///////////////////////////////////////////////
  var areaice = areacalc(elevclassice); //get area of ice-area
  var areasnow = areacalc(elevclassnow); //get area of snow-area
  var areasnowice = areaice.add(areasnow);
  var snowpart = areasnow.divide(areasnowice);
  var voidpart = ee
    .Number(1)
    .subtract(
      ee.Number(areaice).add(ee.Number(areasnow)).divide(areacalc(elevation))
    );

  var ratio1 = areasnow.divide(areaice); //calculate ration between ice and snow
  var ratio2 = areaice.divide(areasnow);
  var multipl = ratio1.multiply(ratio2);

  var sum = ratio1.add(ratio2);
  var ratio3 = ratio1.divide(sum);
  var condition1 = snowpart.gt(0.95); //set condition to 1 (Booliean), when 90% of the snow/ice is covered with Snow
  // if condition TRUE, SnowLine is set to the minimal altitude of the glacier-outline-area (RGI)

  ///////   Get lowest possible SLA   ///////////////////////////////////////////////
  var YYY = function () {
    var lowsla = elevation
      .reduceRegion({
        reducer: ee.Reducer.min(),
        geometry: geometry,
        scale: 30,
        bestEffort: true,
      })
      .get("AVE_DSM");
    return lowsla;
  };

  ///////   Get highest possible SLA   ///////////////////////////////////////////////
  var ZZZ = function () {
    var highsla = elevation
      .reduceRegion({
        reducer: ee.Reducer.max(),
        geometry: geometry,
        scale: 30,
        bestEffort: true,
      })
      .get("AVE_DSM");
    return ee.Number(highsla);
  };

  ///////   Exception Handling   ///////////////////////////////////////////////
  //Set SLA when Ice and Snow Patch do not touch e.o.
  var notouch = function () {
    var lowestsnowelev = function () {
      var numb = elevation
        .clip(maxsnow)
        .reduceRegion({
          reducer: ee.Reducer.min(),
          geometry: geometry,
          scale: 30,
          bestEffort: true,
        })
        .get("AVE_DSM");

      return ee.Number(numb);
    };
    var cond = maxsnow1.size();
    var conditional = ee.Algorithms.If(cond, lowestsnowelev(), ZZZ());
    return ee.Number(conditional);
  };

  var XXX = function (elevclass) {
    var conditionalXXX = ee.Algorithms.If(
      sl1, //check, if SL1 is null or not (e.g. if theres no ice-pixel)
      sl1,
      notouch()
    );
    var sl2 = conditionalXXX;
    return sl2;
  };

  // vectorize the classified map with snow and ice patches
  var classes = forvecfinal.reduceToVectors({
    reducer: ee.Reducer.countEvery(),
    geometry: geometry,
    scale: 20,
    eightConnected: false,
    bestEffort: true,
    maxPixels: 1e9,
  });

  //get the biggest snow and ice patch
  var result = ee.FeatureCollection(classes);
  var filtersnow = ee.Filter.eq("label", 0);
  var maxsnow1 = result.filter(filtersnow).sort("count", false);
  var areasnowvec = ee.Number(maxsnow1.aggregate_sum("count"));
  var maxsnow = ee.Feature(maxsnow1.first());

  var filterice = ee.Filter.greaterThanOrEquals("label", 9);
  var maxice1 = result.filter(filterice).sort("count", false);
  var areaicevec = ee.Number(maxice1.aggregate_sum("count"));
  var maxice2 = ee.Feature(maxice1.first());
  var maxiceifnull = ee.Feature(null).set("count", 0);
  var maxice = ee.Feature(ee.Algorithms.If(maxice2, maxice2, maxiceifnull));

  var fcsnowice = ee.Algorithms.Collection([maxsnow, maxice]);
  var areasnowpatch = ee.Algorithms.If(
    maxsnow1.size(),
    ee.Number(maxsnow.get("count")),
    0
  );
  var areaicepatch = ee.Number(maxice.get("count"));
  var area_all = ee.Number(result.aggregate_sum("count"));

  var areaiceratio = areaicevec.divide(area_all);
  var areasnowratio = areasnowvec.divide(area_all);
  var ratiosnowpatch = ee.Number(areasnowpatch).divide(area_all).multiply(100);
  var ratioicepatch = areaicepatch.divide(area_all).multiply(100);
  var consideredareapatches = ratiosnowpatch.add(ratioicepatch);

  //Extract the zone, where the two patches touches each other
  var img2cl = fcsnowice.reduceToImage(["label"], ee.Reducer.first());
  var biggerice = img2cl.mask(img2cl.select("first").eq(0)).focal_max(2);
  var touchingzone = biggerice.subtract(
    img2cl.mask(img2cl.select("first").gte(9))
  );
  var elevtouch = elevation.addBands(touchingzone, ["first"]);
  var elevationsl = elevtouch.mask(elevtouch.select("first").lt(-1));

  var sl_vector1 = function () {
    var vectorized = touchingzone.reduceToVectors({
      reducer: ee.Reducer.countEvery(),
      geometry: geometry,
      scale: 1,
      maxPixels: 1e8,
      bestEffort: true,
    });
    return vectorized.sort("count", false);
  };

  // calculate mean altitude of the zone where the patches touches
  var meanDictionary = elevationsl.reduceRegion({
    reducer: ee.Reducer.median(),
    maxPixels: 1e8,
    geometry: geometry,
    bestEffort: true,
  });

  var sl1 = meanDictionary.get("AVE_DSM");

  ///////   Exception handling   ///////////////////////////////////////////////
  var conditional1 = ee.Algorithms.If(condition1, YYY(), XXX());
  var sl_vector = ee.Algorithms.If(sl1, sl_vector1(), null);

  var sl = ee.Number(conditional1);
  // calculate std. dev of the extracted Snow line
  var calcstddev = function () {
    var stddev = elevationsl.reduceRegion({
      reducer: ee.Reducer.stdDev(),
      maxPixels: 1e8,
      geometry: geometry,
      bestEffort: true,
    });
    return ee.Number(stddev.get("AVE_DSM"));
  };
  // check if 95% of the glacier is covered with snow --> set sla to the lowest elevation of the glacier
  var check95 = ee.Algorithms.If(condition1, 0, calcstddev());
  var stdDevSLA = ee.Algorithms.If(meanDictionary.size(), check95, 0);

  //add date in MS-Excel readable format
  var imgdate = ee.Number(classified.get("system_time_start"));
  var shortenpre = imgdate.divide(1000).floor(); //the last three digits of the ee.Number are wrong (?)
  var shorten = shortenpre.divide(86400).add(25569); //and for this reason are removed with the division and the floor().

  // Create Feature to return and export informations in a table (e.g. csv)
  var feature = ee.Feature(null);
  feature = feature.set("Snow Cover Ratio", areasnowratio);
  feature = feature.set("SLA MP-Approach", sl);
  feature = feature.set("Ratio Area v/o Snow or Ice", voidpart);
  feature = feature.set("system:time_start", imgdate);
  feature = feature.set("date_MSxlsx", shorten);
  feature = feature.set("otsu", otsu);
  feature = feature.set(
    "Considered Area for MP-SLA Extraction",
    consideredareapatches
  );
  feature = feature.set(
    "stdDev SLA, 0 means no touch between Main patches [m]",
    stdDevSLA
  );
  feature = feature.set("ID_Glacier", glimsid);

  return feature;
};

// hello
