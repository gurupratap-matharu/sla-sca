/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var L4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2"),
  L5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2"),
  L7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2"),
  L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
  L9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"),
  Box =
    /* color: #ffc82d */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
      [
        [
          [-113.31908262896587, 41.80375175369021],
          [-113.31908262896587, 40.6216988523074],
          [-111.80949248004009, 40.6216988523074],
          [-111.80949248004009, 41.80375175369021],
        ],
      ],
      null,
      false
    ),
  DEM = ee.Image("NASA/NASADEM_HGT/001");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
/* Time series visualization using combined Landsat 4-9 C2L2T1 imagery.
   1.  Select and modify the Box to bound selected area
   2.  Run the script
   3.  Modify date range if needed 
   4.  Modify VisParameters (VP) to display best stretch.
   5.  Rerun the script
     
    Final presentation is visualized over a shaded relief from NASA 30m SRTM
    
    --  Script is delivered as is. 
      R. Douglas Ramsey
      Dept. of Wildland Resources
      Quinney College of Natural Resources
      Utah State University #
      Logan, Utah 84322-5230
      doug.ramsey@usu.edu
*/
// ################################################################
// ############### INITIAL PARAMETERS #############################
// Set extent
var fc = Box;
Map.centerObject(fc, 11);
// set start and end years
var startMonth = 7;
var endMonth = 8;
var allYears = ee.List.sequence(1985, 2021); // change years as appropriate.

// Calc location of text
//Casting it to an array makes it possible to slice out the x and y coordinates:
var listCoords = ee.Array.cat(Box.coordinates(), 1);
// get the X-max/min  & Y-max/min coordinates
var xMin = listCoords.slice(1, 0, 1).reduce("min", [0]).get([0, 0]);
var xMax = listCoords.slice(1, 0, 1).reduce("max", [0]).get([0, 0]);
var yMin = listCoords.slice(1, 1, 2).reduce("min", [0]).get([0, 0]);
var yMax = listCoords.slice(1, 1, 2).reduce("max", [0]).get([0, 0]);
// Get distance between two longitudes
var ll = ee.Geometry.Point(xMin, yMin);
var lr = ee.Geometry.Point(xMax, yMin);
var ul = ee.Geometry.Point(xMin, yMax);
var ur = ee.Geometry.Point(xMax, yMax);
var dist = ll.distance(lr).add(ll.distance(ul)).divide(2);
print("dist", dist.divide(100));

//  Calculate point coordinates 0.1 distance from mins
var xPoint = xMin.add(xMax.subtract(xMin).multiply(0.15));
var yPoint = yMin.add(yMax.subtract(yMin).multiply(0.15));

var LabelLocation = ee.Geometry.Point(xPoint, yPoint);

// Visualization parameters.  Modify these as needed.
var VP = {
  opacity: 1,
  bands: ["R", "G", "B"],
  min: [8750, 9463, 8873],
  max: [20000, 16000, 13937],
  gamma: 2,
};
var VP = {
  opacity: 1,
  bands: ["NIR", "R", "G"],
  min: [9735, 7950, 7882],
  max: [19797, 15544, 12904],
  gamma: 2,
}; //st George
var VP = {
  opacity: 1,
  bands: ["R", "G", "B"],
  min: 7355,
  max: 21826,
  gamma: 2,
};
var VP = {
  opacity: 1,
  bands: ["SWIR2", "NIR", "R"],
  min: [7088, 7175, 8091],
  max: [26844, 28366, 30398],
  gamma: 2,
};

// ################################################################
// #################### FUNCTIONS #################################

// TM Filter for clouds
function masksr(image) {
  var qa = image.select("QA_PIXEL");
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = 1 << 4;
  var cloudsBitMask = 1 << 3;
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa
    .bitwiseAnd(cloudShadowBitMask)
    .eq(0)
    .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask);
}

//import the gena package to add annotation to images
var text = require("users/gena/packages:text");

//##########################################################################
//##########################################################################
//################### PROCESS IMAGERY ######################################
//  This sets up the selection and renaming of bands so all images will have the same spectral bands
var l8bands = ["SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7"];
var l57bands = ["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B7"];
var bndNames = ["B", "G", "R", "NIR", "SWIR1", "SWIR2"];

// Get the Image Collections and merge into one
// Load the Landsat 8 ImageCollection.  Mask clouds and add NDVI.
var l9 = L9.filterBounds(fc)
  .map(masksr)
  .map(function (image) {
    return image.select(l8bands).rename(bndNames);
  });
var l8 = L8.filterBounds(fc)
  .map(masksr)
  .map(function (image) {
    return image.select(l8bands).rename(bndNames);
  });
var l7 = L7.filterBounds(fc)
  .map(masksr)
  .map(function (image) {
    return image.select(l57bands).rename(bndNames);
  });
var l5 = L5.filterBounds(fc)
  .map(masksr)
  .map(function (image) {
    return image.select(l57bands).rename(bndNames);
  });
var l4 = L4.filterBounds(fc)
  .map(masksr)
  .map(function (image) {
    return image.select(l57bands).rename(bndNames);
  });

var TMmerged = l4
  .merge(l5)
  .merge(l7)
  .merge(l8)
  .merge(l9)
  .sort("system:time_start");

//  Generate a hillshade to place under the NDVI
var HillShade = ee.Terrain.hillshade(DEM.clip(fc).multiply(2), 135, 45).divide(
  255
);

// set annotation parameters and label location
print("Scale", Map.getScale());
var annotations = [
  {
    position: "left",
    offset: "1%",
    margin: "1%",
    property: "label",
    scale: dist.divide(100),
  }, //Map.getScale() * 2
];

//##########################################################################
//##########################################################################
//  CREATE YEARLY IMAGES, MERGE WITH SHADED RELIEF AND IMBED DATE.
// Create image collection of individual dates. Merge previous year with current
//  year to fill in cloud masked areas if not filled by current year's images.
var Images = ee.ImageCollection.fromImages(
  allYears.map(function (Y) {
    var year = ee.Number(Y);
    var image1 = TMmerged.filter(
      ee.Filter.calendarRange(year.subtract(1), year.subtract(1), "year")
    )
      .filter(ee.Filter.calendarRange(startMonth, endMonth, "month"))
      .median()
      .clip(fc);
    var image2 = TMmerged.filter(ee.Filter.calendarRange(Y, Y, "year"))
      .filter(ee.Filter.calendarRange(startMonth, endMonth, "month"))
      .median()
      .set("year", Y)
      .clip(fc);

    return image1.blend(image2).set("year", Y);
  })
);
print(Images);

//Blend the Hillshade, NDVI,water, and  state boundaries for every image in the collection.
var sequence = Images.map(function (image) {
  var year = image.get("year");
  var theImage = image
    .visualize(VP)
    .multiply(HillShade)
    .set({ label: ee.Number(year).int() });
  var annotated = text
    .annotateImage(theImage, {}, LabelLocation, annotations)
    .set("year", year);
  return annotated;
});

// Display a selected year.
var image = sequence.filter(ee.Filter.eq("year", 2021)); //.first()
Map.addLayer(Images.filter(ee.Filter.eq("year", 2021)), VP, "Landsat");
//Map.addLayer(image,{},'Image')

// Define animation arguments.
var thumbArgs = {
  dimensions: 1500,
  region: fc,
  crs: "EPSG:3857", //'EPSG:3857',26912
};

// Display the thumbnail.
print(ui.Thumbnail(image, thumbArgs));

// ################################################################
// ################## ANIMATED THUMBNAIL ##########################
// ################################################################

//Create an animated thumbnail
// Visualization parameters.
var args = {
  crs: "EPSG:3857",
  dimensions: "800",
  region: fc,
  framesPerSecond: 4,
};

// Create the animated thumbnail and add it to the map.
var thumb = ui.Thumbnail({
  image: sequence,
  params: args,
  style: {
    position: "top-right",
    width: "800 px",
  },
});
Map.add(thumb);

// hello
