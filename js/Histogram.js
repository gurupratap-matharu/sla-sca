/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var RGI = ee.FeatureCollection("users/josiaszeller/RGVI_v6"),
  alos = ee.Image("JAXA/ALOS/AW3D30/V2_2"),
  srtm = ee.Image("USGS/SRTMGL1_003");
/***** End of imports. If edited, may not auto-convert in the playground. *****/

exports.sla_extract = function (img) {
  ///////   Creating variables   ///////////////////////////////////////////////
  var glimsid = img.get("GLIMSID");
  var deminfo = img.get("deminfo");
  var geometry = RGI.filterMetadata("GLIMSId", "equals", glimsid);
  var classified = img;

  // load DEM of the Glacier to extract the height distribution of the glacier

  var dem1 = srtm.select("elevation").rename("AVE_DSM").clip(geometry); //rename srtm dem to match the alos dem, clip to geometry
  var dem2 = alos.select("AVE_DSM").clip(geometry);
  var demselector = ee.Algorithms.IsEqual(
    ee.String(deminfo),
    ee.String("ALOS")
  ); //boolean to check, if dem selection == ALOS
  var demselection = ee.Image(ee.Algorithms.If(demselector, dem2, dem1)); //select the DEm that should be used
  var demglacier = demselection; //assignment

  //////////////    CALCULATE AREA OF IMAGE AND RETURN AREA   ////////////////
  var areacalc = function (image4) {
    //function to calculate the area of an image (area of a band without masked pixels)
    var prep = image4.select("AVE_DSM"); //take one bands (here green) to calculate the area.
    var pixelarea = prep.multiply(ee.Image.pixelArea());
    var areacalc = pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 30,
      maxPixels: 1e10,
    });
    var numb = areacalc.getNumber("AVE_DSM");
    return numb; //return the number of the area
  };
  //////////////    START CODE    ///////////////////////////////////////////////
  var dcmask = classified.select("classification").neq(8);
  var mask1 = classified.select("classification").neq(7);
  var mask2 = classified.select("classification").neq(4);
  var mask3 = classified.select("classification").neq(3);
  var mask4 = classified.select("classification").neq(2);
  var masksnow = classified.select("classification").neq(6);
  var classmask1 = classified.mask(dcmask).clip(geometry);
  var classmask2 = classmask1.updateMask(mask1).clip(geometry);
  var classmask3 = classmask2.updateMask(mask2).clip(geometry);
  var classmask4 = classmask3.updateMask(mask3).clip(geometry);
  var classmask5 = classmask4.updateMask(mask4).clip(geometry);
  var classmask6 = classmask5.updateMask(masksnow).clip(geometry); //.unmask(10).clip(geometry)

  var snowall = classified.mask(classified.select("classification").eq(6)); //combine snow with shadow on snow
  var snowshadowlayer2 = snowall.set("classification", 1);

  var snowunited = ee.ImageCollection([classmask6, snowshadowlayer2]);
  var twoclasses = snowunited.reduce(ee.Reducer.max());

  ///////////////////////////////////////////// ELEVATION ANALYSIS ////////////////////////
  var elevation = demglacier.select("AVE_DSM").clip(geometry);

  var lowsla1 = elevation
    .reduceRegion({
      reducer: ee.Reducer.min(),
      geometry: geometry,
      scale: 30,
    })
    .get("AVE_DSM");

  var highsla = elevation
    .reduceRegion({
      reducer: ee.Reducer.max(),
      geometry: geometry,
      scale: 30,
    })
    .get("AVE_DSM");

  ///////   Calculations to export SCR   ///////////////////////////////////////////////
  var elevclassice = elevation.mask(
    classified.select("classification").clip(geometry).eq(0)
  );
  var elevclassnow = elevation.mask(
    classified.select("classification").clip(geometry).eq(1)
  );

  var areaice = areacalc(elevclassice);
  var areasnow = areacalc(elevclassnow);

  var sla_normal = function () {
    var areasnow = areacalc(elevclassnow);
    var areasnowice = areaice.add(areasnow);
    var snowpart = areasnow.divide(areasnowice);

    var ratio1 = areasnow.divide(areaice);
    var ratio2 = areaice.divide(areasnow);
    var multipl = ratio1.multiply(ratio2);

    var sum = ratio1.add(ratio2);
    var ratio3 = ratio1.divide(sum);
    var condition1 = snowpart.gt(0.9);
    var condition2 = multipl.eq(1);

    var YYY = function () {
      var lowsla = elevation
        .reduceRegion({
          reducer: ee.Reducer.min(),
          geometry: geometry,
          scale: 30,
        })
        .get("AVE_DSM");
      return lowsla;
    };

    var elevclass = elevclassice.addBands(elevclassnow);

    ///////   Creat histogram for Ice and SNow   ///////////////////////////////////////////////
    var XXX = function (elevclass) {
      var hist_snow = elevation
        .updateMask(classified.select("classification").eq(1))
        .reduceRegion({
          //calculate histogram for snow
          reducer: ee.Reducer.autoHistogram(40, 50),
          scale: 30,
          geometry: geometry,
          maxPixels: 1e13,
        });

      var hist_snow_array = ee.Array(hist_snow.get("AVE_DSM")).toList(); //convert the Object (list) to a list with features
      var features_snow = hist_snow_array.map(function (item) {
        item = ee.List(item);
        return ee.Feature(null, {
          elevation: item.get(0),
          countsnow: item.get(1),
        });
      });

      var hist_ice = elevation
        .updateMask(classified.select("classification").eq(0))
        .reduceRegion({
          //calculate histogram for ice. classification 0 is ice
          reducer: ee.Reducer.autoHistogram(40, 50),
          geometry: geometry,
          scale: 30,
          maxPixels: 1e13,
        });

      var hist_ice_array = ee.Array(hist_ice.get("AVE_DSM")).toList(); //convert the Object (list) to a list with features
      var features_ice = hist_ice_array.map(function (item) {
        item = ee.List(item);
        return ee.Feature(null, {
          elevation: item.get(0),
          countice: item.get(1),
        });
      });

      var fcFilter = ee.Filter.equals({
        // create filter to use it for the inner join of the two FeatureCollections (snow and ice)
        leftField: "elevation",
        rightField: "elevation",
      });

      var fcFilter1 = ee.Filter.equals("elevation"); //join both FeatureCollections
      var joiner = ee.Join.inner();
      var fcJoin = joiner.apply(features_snow, features_ice, fcFilter);

      ///////   Calculate ratio to extract SLA   ///////////////////////////////////////////////
      var calcfc = fcJoin.map(function (i) {
        //map function to calculate the ration between sum of pixels and the difference of pixels over featurecollection
        var snowfeat = ee.Feature(i.get("primary"));
        var icefeat = ee.Feature(i.get("secondary"));
        var snow = snowfeat.getNumber("countsnow");
        var ice = icefeat.getNumber("countice");
        var dif = snow.subtract(ice).abs(); //calculate the difference per elevation
        var sum = snow.add(ice); //sum of ice and snow pixels per elevation
        var ratio = sum.divide(dif); //ratio between sum and difference per elevation stept
        var snowfeatdif = snowfeat
          .set("dif", dif)
          .set("countice", ice)
          .set("ratio", ratio)
          .set("sum", sum); //write all properties in one feature
        return ee.Feature(snowfeatdif); //return feature and create featurecollection with one element per elevation step
      });

      var lengthnew = ee.Number(calcfc.size().multiply(0.4).floor());
      var meansum = ee.Number(calcfc.aggregate_mean("sum"));
      var meanice = ee.Number(calcfc.aggregate_mean("ice"));

      var filteredfc = calcfc.limit(lengthnew, "sum", false);
      var sorted = filteredfc.sort("ratio", false); //set to true if you want the minimum in difference where the sum of snow and ice pixels per elevation is maximised

      var sl = ee.Number(sorted.first().get("elevation"));
      return sl;
    };

    var conditional1 = ee.Number(ee.Algorithms.If(condition1, YYY(), XXX()));

    return conditional1;
  };

  var calcarea1 = function (image4) {
    //function to calculate the area of an image (area of a band without masked pixels)
    var prep = image4.select("classification"); //take one bands (here green) to calculate the area.
    var pixelarea = prep.multiply(ee.Image.pixelArea());
    var areacalc = pixelarea.reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: 30,
      maxPixels: 1e10,
    });
    var numb = areacalc.getNumber("classification");
    return ee.Number(numb); //return the number of the area
  };

  var ifnosnow = function () {
    var condition = calcarea1(snowshadowlayer2);
    var sla = ee.Algorithms.If(condition, sla_normal(), ee.Number(highsla));
    return ee.Number(sla);
  };

  var condslafin = calcarea1(
    classified.select("classification").clip(geometry).eq(0)
  );

  var sla_final = ee.Algorithms.If(areaice, ifnosnow(), ee.Number(lowsla1)); //check, if theres no ice-pixel. If no ice is detectes, lowSLA is returned
  var sl = sla_final;

  ///////   Date in MSExcel Readable Format   ///////////////////////////////////////////////
  var imgdate = ee.Number(classified.get("system_time_start"));
  //add date in a Microsoft Excel usable format
  var shortenpre = imgdate.divide(1000).floor(); //the last three digits of the ee.Number are wrong (?)
  var shorten = shortenpre.divide(86400).add(25569); //and for this reason are removed with the division and the floor().

  var feature = ee.Feature(null);
  feature = feature.set("system:time_start", imgdate);
  feature = feature.set("SLA_Histogram_Approach", sl);
  feature = feature.set("date_MSxlsx", shorten);

  return feature;
};

// Hello
