/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var ING = ee.FeatureCollection(
  "users/lcsruiz/Mapping_seasonal_glacier_melt_across_the_ANDES_with_SAR/Glaciares_Arg_Andes_dissolve"
);
/***** End of imports. If edited, may not auto-convert in the playground. *****/

//CloudScore originally written by Matt Hancher and adapted for S2 data by Ian Housman
//Adapted for Landsat 5 to 8 and Sentinel 2 data by Josias Zeller

exports.CloudShadows = function (img) {
  ///////   Set Parameters   ///////////////////////////////////////////////
  var cloudThresh = 20; //Ranges from 1-100.Lower value will mask more pixels out. Generally 10-30 works well with 20 being used most commonly
  var cloudHeights = ee.List.sequence(0, 500, 25); //Height of clouds to use to project cloud shadows
  var irSumThresh = 0.45; //org 0.35Sum of IR bands to include as shadows within TDOM and the shadow shift method (lower number masks out less)
  var dilatePixels = 15; //Number of pixels to buffer cloud and cloud shadow mask by
  var contractPixels = 10; //Number reduce cloud and cloud shadow mask by to reduce inclusion of single-pixel comission errors

  var ing_id = img.get("ING_ID");
  var geometry = ING.filterMetadata("ID_local", "equals", ing_id);

  //////////////////////////////////////////////////////////////////////////
  // FUNCTIONS
  //////////////////////////////////////////////////////////////////////////
  var rescale = function (img, exp, thresholds) {
    return img
      .expression(exp, { img: img })
      .subtract(thresholds[0])
      .divide(thresholds[1] - thresholds[0]);
  };

  ////////////////////////////////////////////////////////////////////////////////
  // Cloud masking algorithm for Sentinel2
  //Built on ideas from Landsat cloudScore algorithm
  //Currently in beta and may need tweaking for individual study areas
  function CloudScore(img) {
    // Compute several indicators of cloudyness and take the minimum of them.
    var score = ee.Image(1);

    // Clouds are reasonably bright in the blue band.
    score = score.min(rescale(img, "img.blue", [0.1, 0.5]));

    // Clouds are reasonably bright in all visible bands.
    score = score.min(
      rescale(img, "img.red + img.green + img.blue", [0.2, 0.8])
    );

    //Clouds are moist
    var ndmi = img.normalizedDifference(["nir", "swir1"]);
    score = score.min(rescale(ndmi, "img", [-0.1, 0.1]));

    // However, clouds are not snow.
    var ndsi = img.normalizedDifference(["green", "swir1"]);
    score = score.min(rescale(ndsi, "img", [0.4, 0.1]));

    score = score.multiply(100).byte();

    return img.addBands(score.rename("cloudScore"));
  }

  //////////////////////////////////////////////////////////////////////////
  /*
   * Implementation of Basic cloud shadow shift
   *
   * Author: Gennadii Donchyts
   * License: Apache 2.0
   */
  function projectShadows(cloudMask, image, cloudHeights, dilatePixels) {
    var meanAzimuth = image.get("SUN_AZIMUTH");
    var meanZenith1 = image.get("SUN_ELEVATION");
    var rightangle = ee.Number(90);
    var meanZenith = rightangle.subtract(meanZenith1);
    ///////////////////////////////////////////////////////

    //Find dark pixels
    var darkPixels = image
      .select(["nir", "swir1", "swir2"])
      .reduce(ee.Reducer.sum())
      .lt(irSumThresh); //.gte(1);

    //Get scale of image
    var nominalScale = cloudMask.projection().nominalScale();

    //Find where cloud shadows should be based on solar geometry
    //Convert to solar geometry to radians
    var azR = ee
      .Number(meanAzimuth)
      .multiply(Math.PI)
      .divide(180.0)
      .add(ee.Number(0.5).multiply(Math.PI));
    var zenR = ee
      .Number(0.5)
      .multiply(Math.PI)
      .subtract(ee.Number(meanZenith).multiply(Math.PI).divide(180.0));

    //Find the shadows
    var shadows = cloudHeights.map(function (cloudHeight) {
      cloudHeight = ee.Number(cloudHeight);

      var shadowCastedDistance = zenR.tan().multiply(cloudHeight); //Distance shadow is cast
      var x = azR
        .cos()
        .multiply(shadowCastedDistance)
        .divide(nominalScale)
        .round(); //X distance of shadow
      var y = azR
        .sin()
        .multiply(shadowCastedDistance)
        .divide(nominalScale)
        .round(); //Y distance of shadow
      return cloudMask.changeProj(
        cloudMask.projection(),
        cloudMask.projection().translate(x, y)
      );
    });

    var shadowMask = ee.ImageCollection.fromImages(shadows).max();

    //Create shadow mask
    shadowMask = shadowMask.and(cloudMask.not());
    shadowMask = shadowMask.focal_min(contractPixels).focal_max(dilatePixels);

    shadowMask = shadowMask.and(darkPixels);

    var cloudShadowMask = shadowMask; //.or(cloudMask);

    image = image.addBands(shadowMask.rename(["cloudShadow"])); //.updateMask(cloudShadowMask.not()).addBands(shadowMask.rename(['cloudShadow']));

    return image;
  }

  //////////////////////////////////////////////////////////////////////////
  // ACTIONS
  ////////////////////////////////////////////////////////////////////////////////
  // running the cloud score
  img = CloudScore(img);
  var cloudScore = img.select("cloudScore");

  ////////////////////////////////////////////////////////////////////////////////
  // masking the cloud score
  var CloudMask = cloudScore
    .gt(cloudThresh)
    .focal_min(contractPixels)
    .focal_max(dilatePixels);
  var CloudMasked = img.updateMask(CloudMask.not());

  ////////////////////////////////////////////////////////////////////////////////
  // projecting the cloud shadow
  var CloudShadowMasked = projectShadows(
    CloudMask,
    img,
    cloudHeights,
    dilatePixels
  );
  var Cloudshadowclipped = CloudShadowMasked.clip(geometry);

  return Cloudshadowclipped;
};

// Hello
