/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var alos = ee.Image("JAXA/ALOS/AW3D30/V2_2"),
    srtm = ee.Image("USGS/SRTMGL1_003"),
    ING = ee.FeatureCollection("users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
exports.ABextraction = function(imgprep) {  

// script original de Zeller modificado para usar el ING como base para el analisis

// Value in meters indicates the height interval
// If the height is set too high, it is possible that the SLA in the upper section will no longer work.
// The current function can only cover a SLA up to the top three height levels.

///////   Creating variables   ///////////////////////////////////////////////
var img = imgprep;
var otsu = img.get('otsu');
var bin = ee.Number(25); // bin in Meters
var ing_id =     img.get('ING_ID');
var deminfo =     img.get('deminfo');
var geometry =    ING.filterMetadata('ID_local','equals',ing_id);



// load DEM of the Glacier to extract the height distribution of the glacier

var dem1 = srtm.select('elevation').rename('AVE_DSM').clip(geometry); //rename srtm dem to match the alos dem, clip to geometry
var dem2 = alos.select('AVE_DSM').clip(geometry);
var demselector = ee.Algorithms.IsEqual(ee.String(deminfo), ee.String('ALOS')); //boolean to check, if dem selection == ALOS
var demselection = ee.Image(ee.Algorithms.If(demselector, dem2, dem1)); //select the DEm that should be used
var demglacier = demselection; //assignment

// store the lowest point of the glacier (height)
var lowestpoint = ee.Number(demglacier.reduceRegion({
      reducer: ee.Reducer.min(),
      geometry: geometry,
      scale : 30,
      bestEffort: true
      }).get('AVE_DSM'));
//store the highest point of the glacier
var highestpoint = ee.Number(demglacier.reduceRegion({
      reducer: ee.Reducer.max(),
      geometry: geometry,
      scale : 30,
      bestEffort: true
      }).get('AVE_DSM'));
///////   Start Code   ///////////////////////////////////////////////

//create list with the levels (bin height can be set in the parameter part)
var altitudelevel =   ee.List.sequence(lowestpoint,highestpoint, bin);
//Add DEM as a Band to the classified image
var dummy   =         ee.Image(1).select(['constant'], ['area']);
var image_dem1 =       img.addBands(demglacier);
var image_dem =       image_dem1.addBands(dummy).clip(geometry);   


//function to calculate the area of an image (area of a band without masked pixels) 
var areacalc = function(image4) { 
                      var prep = image4.select('area'); //take one bands (here AVE) to calculate the area. 
                     //var pixelarea = prep.multiply(ee.Image.pixelArea());
                      var areacalc1 = prep.reduceRegion ({
                                             reducer: ee.Reducer.sum(),
                                             geometry: geometry,
                                             scale: 30,
                                             maxPixels: 1e11,
                                             bestEffort: true
                                        }) ;
                      var numb = areacalc1.getNumber('area') ;
                      return numb; //return the number of the area
                      };

// function to get the ration of Snow/Ice distribution on the glacier for each bin
//input is the lowest height (masl) of the bin
// void is everything but snow and ice (water, shadow, debris cover, ...)
var getlevelinfo = function(level) {
  var binpre =      bin.subtract(0.001);                   
  var lowerb =      ee.Number(level);
  var higherb =     ee.Number(level).add(binpre);
// filter the image with the bin boundaries
  var lowermask =       image_dem.select('AVE_DSM').lte(higherb);
  var highermask =      image_dem.select('AVE_DSM').gte(lowerb);
// mask the image with the created mask  
///--------------
var areacalc1 = function(image4) { 
                      var prep = image4.select('area'); //take one bands (here AVE) to calculate the area. 
                    // var pixelarea = prep.multiply(ee.Image.pixelArea());
                      var areacalc2 = prep.reduceRegion ({
                                             reducer: ee.Reducer.sum(),
                                             geometry: geometry,
                                             scale: 30,
                                             maxPixels: 1e11,
                                             bestEffort: true
                                        }) ;
                      var numb = areacalc2.getNumber('area') ;
                      return numb; //return the number of the area
                      };

  var classedlevel1 =   image_dem.mask(lowermask).updateMask(highermask);

  var filter_masked =   classedlevel1.select('AVE_DSM').gt(0);
  var filter_snow =     classedlevel1.select('classification').eq(1); //select snow
  var filter_ice =      classedlevel1.select('classification').eq(0); //select ice
  var filter_nothing =      classedlevel1.select('classification').neq(11); //select ice
  var classedlevel =    classedlevel1.mask(filter_masked)//.clip(geometry); //.select('classification').add(1);
  var snowlevel =       classedlevel.mask(filter_snow)//.clip(geometry); //.select('classification').add(1);
  var icelevel =        classedlevel.mask(filter_ice)// .clip(geometry); //.select('classification').add(1); //add 1, that the AreaCalculation works
  var calcvoidarea =    classedlevel.mask(filter_nothing)//.clip(geometry);
  var arealevel =       ee.Number(areacalc1(calcvoidarea)); //calculate area of the whole bin
  var areasnow =        ee.Number(areacalc1(snowlevel)); //calculate area of the snow cover
  var areaice =         ee.Number(areacalc1(icelevel)); //calculate area of the ice cover
 
  var areavoid =        arealevel.subtract(areasnow).subtract(areaice); //calculate area of the void cover
 // calculate the rations to distribute the void-part in the same proportion as snow and ice
 // in the corresponding bin
  var snowratio =       areasnow.divide(arealevel);
  var iceratio =        areaice.divide(arealevel);
  var voidratio =       ee.Number(1).subtract(snowratio).subtract(iceratio);

  var totalsnowice =      areaice.add(areasnow);
  
  var totalice =          areaice.add((areaice.divide(totalsnowice)).multiply(areavoid));
  var totalsnow =         areasnow.add((areasnow.divide(totalsnowice)).multiply(areavoid));
  
  var totarea =        totalice.add(totalsnow);
  var ratiototalsnow =      (areasnow.add((areasnow.divide(totalsnowice)).multiply(areavoid))).divide(arealevel);

///////   Create Feature   ///////////////////////////////////////////////
  var feature = ee.Feature(null);
feature = feature.set('arealevel_total', totarea);
feature = feature.set('areavoid', areavoid); 
feature = feature.set('snowtotal',  totalsnow);
feature = feature.set('ratio totalsnow (SCR)',  ratiototalsnow);
feature = feature.set('level lower boundary',   level);
return(feature);
};

var FC =  ee.FeatureCollection(altitudelevel.map(getlevelinfo));

///////   SLA und SCR Calculations   ///////////////////////////////////////////////
var totalareaglacier =  ee.Number(areacalc(image_dem));
var totalareasnow =     FC.aggregate_sum('snowtotal');
var totalareavoid =     FC.aggregate_sum('areavoid');

var ratiosnow =         totalareasnow.divide(totalareaglacier);
var voidarearatio =     totalareavoid.divide(totalareaglacier);

///////   Get Time to add to the Output   ///////////////////////////////////////////////
var syststart =   ee.Number(img.get('system_time_start'));
//add date in a Microsoft Excel usable format 
var shortenpre =  syststart.divide(1000).floor(); //the last three digits of the ee.Number are wrong (?)
var shorten =     shortenpre.divide(86400).add(25569); //and for this reason are removed with the division and the floor().

///////   Extract SLA   ///////////////////////////////////////////////

var countlist = (FC.size()); //get Number of Features in the FeatureCollection
// Creat list from the lowest Point with a Step of the level height
// the subtraction of four (4) is to shorten the list by the last four levels
var list =          ee.List.sequence(lowestpoint,null,bin,(countlist.subtract(4)));

var addsla =    function(nmbr){
var startnmr =      ee.Number(nmbr);
var level =         startnmr;
var first =         ee.Number((FC.filterMetadata('level lower boundary','equals', startnmr).first()
                            .get('ratio totalsnow (SCR)')))
                            .gt(0.5);
var secondfilter =  startnmr.add(bin);
var thirdfilter =   secondfilter.add(bin);
var second =        ee.Number((FC.filterMetadata('level lower boundary','equals', secondfilter).first()
                            .get('ratio totalsnow (SCR)')))
                            .gt(0.5);
var third =         ee.Number((FC.filterMetadata('level lower boundary','equals', thirdfilter).first()
                            .get('ratio totalsnow (SCR)')))
                            .gt(0.5);
var sumcheck =      first.add(second).add(third).eq(3);
var feature =       ee.Feature(null)
                            .set('level', nmbr)
                            .set('check', sumcheck)
                            .set('level lower boundary', level);
return              feature;
                    };
                    

//Map the function (addsla) over the list
var slafcpre = list.map(addsla);
//filter the FeatureCollection to get a list with all feature that fulfill the condition
// that the three levels above are also covered with at least 50% snow
var slafc = slafcpre.filter(ee.Filter.eq('check',1));
//select the first Feature to extract the SLA
var slaf = ee.Feature(slafc.get(0));

//check if there's is a feature to reject the inquiry if there's no feature
var condition = slafc.size();

// return the lower level boundary from the feature or return -9999 if the condition is not fulfilled
var slacondition = ee.Number(ee.Algorithms.If(condition, ee.Number(slaf.get('level lower boundary')), ee.Number(-9999)));

//check if no bin is covered bis <50% of snow - set SLA to the lowest point on the glacier


///////   Exception handling   ///////////////////////////////////////////////
var alternativesla = function() {
var altsla01 = (FC.filterMetadata('ratio totalsnow (SCR)', "greater_than", 0.4999999 )).size();
var altsla02 = FC.filterMetadata('ratio totalsnow (SCR)', "greater_than", 0.4999999)
                  .sort('ratio totalsnow (SCR)',true)
                  .first()
                  .get('level lower boundary');
var altsla = ee.Number(ee.Algorithms.If(altsla01,altsla02,lowestpoint));
return altsla;
};

var slafeat = function(){
var sla = ee.Number(slaf.get('level lower boundary'));
return sla;
};


///////   Set SLA   ///////////////////////////////////////////////
var sla = slafeat(); 


/////////////////////////
var feature = ee.Feature(null);
feature = feature.set('Snow Cover Ratio (Approximation)',      ratiosnow);
feature = feature.set('SLA AB-Approach',                    sla );
feature = feature.set('Ratio Area w/o Snow or Ice',  voidarearatio);
feature = feature.set('system:time_start',      syststart);
feature = feature.set('date_MSxlsx',            shorten);
feature = feature.set('otsu',                   otsu);
feature = feature.set('ID_Glacier',             ing_id);
var feattoexp = feature;

return feattoexp;
};

// Hello