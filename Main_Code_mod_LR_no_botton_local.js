/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var RGI = ee.FeatureCollection("users/josiaszeller/RGVI_v6"),
    ING = ee.FeatureCollection("users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve");
/***** End of imports. If edited, may not auto-convert in the playground. *****/



// Version 3 - Beta
// Questions? -> info@josiaszeller.ch
// 2020-11-29

// Modificado por Lucas Ruiz

//////////////////     NEW FEATURES ///////////////////////////

// DEM-Selection SRTM and ALOS)
// Improved Legend in MapView
// GeoTiff-Download for the classification can be activated on Line 175

///////////////////////////////////////////////////////////////
/*
Map.addLayer(RGI,{color: 'green'}, 'Randolph Glacier Inventory 6.0');
Map.onClick(function(coords) {  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  print('GLIMS ID: ',ee.String(RGI.filterBounds(click_point).first().get('GLIMSId')));
}) ;

*/
Map.addLayer(ING,{color: 'green'}, 'ING (2018)');
Map.onClick(function(coords) {
  var location = 'lon: ' + coords.lon.toFixed(4) + ' ' +
                 'lat: ' + coords.lat.toFixed(4);
  var click_point = ee.Geometry.Point(coords.lon, coords.lat);
  print('ING ID: ',ee.String(ING.filterBounds(click_point).first().get('ID_local')));
}) ;

//var getglims = ui.Textbox('Enter GLIMS ID here');
//print(getglims);
//print('Not sure what the GLIMS ID is?');
//print('Just click on a glacier');


// Boton RUN
//var runbutton = ui.Button({
//    label: "Run",
//    onClick: function() {

// Print Parameterization to the Console
//var startyearslider =   ui.Slider(1984,2021, 2019,1);
// endyearslider =     ui.Slider(1984,2021, 2019,1);
//var startdoy =          ui.Slider(1,365, 200,1); 
//var enddoy =            ui.Slider(1,365, 210,1);
//var cloudc =            ui.Slider(0,100, 80);
//var hsstartslider =     ui.Slider(0, 365, 1);
//var hsendslider =       ui.Slider(0, 365, 2);
//var demselect =         ui.Select(['ALOS', 'SRTM']);
//print('Select an evaluation period');
//print(startyearslider,  'Year Start');
//print(endyearslider,    'Year End');
//print(startdoy,         'Start of Period (DOY)');
//print(enddoy,           'End of Period (DOY)');
//print(cloudc,           'Cloud cover [%]');
//print(hsstartslider,    'Start <add Hill Shadow> (DOY)');
//print(hsendslider,      'End <add Hill Shadow> (DOY)');
//print(demselect,        'Select a DEM');


// boton Crear image collection
//var button = ui.Button({
//  label: 'Create Image Collection',
//  onClick: function() {

//var glimsid = getglims.getValue();
//var startyear =           startyearslider.getValue();
//ar endyear =             endyearslider.getValue();
//var filterDOYstart =      startdoy.getValue();    // Start of Range: use DOY-Number to filter the range
//var filterDOYend =        enddoy.getValue();      // End of Range: DOY-Number to filter the range
//var cloudiness =          cloudc.getValue();      //% (Area over the Glacier that is not covered with Clouds)
//var coverage =            90; //percentage of the glacier area, that has to be covered by the image
//var hsstart =             ee.Number(hsstartslider.getValue());
//ar hsend =               ee.Number(hsendslider.getValue());
//var dem =                 demselect.getValue();


// Modificado por Lucas Ruiz
// Parametros de entrada seteados para glaciar Alerce
var glimsid = 'G288174E41166S';          // PERITO MORENO = 'G286783E50561S';
var startyear =           ee.Number(2021);
var endyear =             ee.Number(2021);
var filterDOYstart =      ee.Number(90);    // Start of Range: use DOY-Number to filter the range
var filterDOYend =        ee.Number(100);      // End of Range: DOY-Number to filter the range
var cloudiness =          ee.Number(80);// cloudc.getValue();      //% (Area over the Glacier that is not covered with Clouds)
var coverage =            90; //percentage of the glacier area, that has to be covered by the image
var hsstart =             ee.Number(1);
var hsend =               ee.Number(2);
var dem =  'SRTM';//               demselect.getValue();

// Hello

// MODIFICADO LUCAS RUIZ 
// imprimir en consola las variables seleccionadas
print(glimsid,'RGI GLACIAR ALERCE')
print(startyear,  'Year Start');
print(endyear,    'Year End');
print(filterDOYstart,         'Start of Period (DOY)');
print(filterDOYend,           'End of Period (DOY)');
print(cloudiness,           'Cloud cover [%]');
print(coverage,           'percentage of the glacier area coverage by the image [%]');
print(hsstart ,    'Start <add Hill Shadow> (DOY)');
print(hsend,      'End <add Hill Shadow> (DOY)');
print(dem,        'Selected DEM');











///////////////////   VISUALISATION   ////////////////////////////////////



var palette = ['f27324', '24a3f2', '00ff00', 'ff00cc','7e23b5','ffffff',  'c39797','000000', 'e5ffe5'];
var geometry =            RGI.filterMetadata('GLIMSId','equals',glimsid);
                          Map.centerObject(geometry);


print(geometry);
print(RGI.first());



///////////////////   MODULE PREP    /////////////////////////////////////
var composites =          require('users/lcsruiz/OPTICOS:Preparation_Module');
var add_Cloudshadow =     require('users/lcsruiz/OPTICOS:Cloudshadow_Mask');

var hsboolean =           0;
var resultPRE =           composites.imageCollection(glimsid, startyear, endyear, cloudiness, coverage, filterDOYstart,filterDOYend, hsboolean, dem);
var result =              resultPRE.filterMetadata('system:time_start','not_equals', 1443176819706); //filter out this image due to wrong scaled bands
var prepfinal1 =          result.map(add_Cloudshadow.CloudShadows);

///////////////////   Add Hill Shadow Information   //////////////////////
//Split the Collection into Parts to write a 0 or 1
var filterHSY =         ee.Filter.calendarRange(hsstart, hsend,'day_of_year');
var filterHSN =         ee.Filter.calendarRange(hsend, hsstart, 'day_of_year');
// add a Yes to the 'hsboolean' property (hill shadow boolean)
var imgcolY =   prepfinal1.filter(filterHSY).map(function(img)
                {var img1 = img.set('hsboolean',1);
                  return img1;
                });
// add a No to the 'hsboolean' property (hill shadow boolean)
var imgcolN =   prepfinal1.filter(filterHSN).map(function(img)
                {var img1 = img.set('hsboolean',0);
                  return img1;
                });
// Merge the to parts and sort them (ascending, starting with the start of the year.
var prepfinal = ee.ImageCollection(imgcolN.merge(imgcolY).sort('system:time_start'));
var imagecount = prepfinal.size();
//print (imagecount)
print('Number of available images:', imagecount);
print('Calculation stops after 5 minutes!');
    //print(prepfinal);
Map.addLayer(prepfinal.first().reproject((prepfinal.first().select('blue').projection())), {bands: ['red','green', 'blue'], min:0, max:1},'First available Image in this Collection');

print('Classification can be started directly');

// boton start classification
//var buttonDT = ui.Button({
//  label: 'Start Classification',
//  onClick: function() {
print('Switch to the Task-Tab to Export Data to Drive');
//comment out this classifier and activate the RF-Classifier
var decisiontree = require('users/lcsruiz/OPTICOS:Classifier');


var mapCollection = prepfinal.map(decisiontree.decisiontree);
// set position of panel
var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
// Create legend title
var legendTitle = ui.Label({
  value: 'Legend',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette = ['f27324', '24a3f2', '00ff00', 'ff00cc','7e23b5',  'c39797', 'e5ffe5', '00ffff', '00cccc' ];
var palette1 = ['f27324', '24a3f2', '00ff00', 'ff00cc','7e23b5',           'e5ffe5',            '00cccc'];
//                0         1         2         3        4           5       6          7         8
// name of the legend
var names = ['Ice', 'Snow', 'Water', 'Debris cover','Clouds',  'Shadow on Snow', 'Undefined Shadow'];

// Add color and and names
for (var i = 0; i < 6; i++) {
  legend.add(makeRow(palette1[i], names[i]));
  }  
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);
Map.addLayer(mapCollection.first('false').select('classification'), {min:0, max:8, palette:palette}, 'Classification of the last Image in the Collection');


// -------------------------------------------------------------------------
// Export the Collection MODIFICADO LUCAS RUIZ
// -------------------------------------------------------------------------
// aqui va codigo para expotar la coleccion de imagenes clasificadas como tiff
// Por ahora solo exporta la ultima imagen de la coleccion

var exampletif = mapCollection.first().select('classification');

Export.image.toDrive({image: exampletif, scale: 20, description: 'imageExport_Alerce19'});


print(mapCollection,'mapCollection')
// termina export Lucas Ruiz



var getAB =               require('users/lcsruiz/OPTICOS:Altitude_Bins');
var mp =                  require('users/lcsruiz/OPTICOS:Main_Patches');
var histogram =           require('users/lcsruiz/OPTICOS:Histogram');

var ABcollection =                mapCollection.map(getAB.ABextraction);
Export.table.toDrive(ABcollection, 'AltitudeLevelsApproach_ALERCE19', 'slamontitoringdata_test_Alerce');
print(ABcollection,'Altitude level collection')
//var Histoocollection =               mapCollection.map(histogram.sla_extract);
//Export.table.toDrive(Histoocollection, 'Histogram', 'testing1');

var MPcollection =                mapCollection.map(mp.sla_extract_patch);
Export.table.toDrive(MPcollection,'MainPatchesApproach_ALERCE19', 'slamontitoringdata_test_Alerce');
print(MPcollection,'Main Patch collection')
//print(MPcollection)
//}});
//print(buttonDT);
//}

//});
//print(button);

// END Boton RUN
//}});
//print(runbutton);

// *   Every module can be accessed: Just add 'https://code.earthengine.google.com/?scriptPath='
// *   to the link in this Script (string that is called with the <require>-command)

//-----------------------------------------------------------------------------------
