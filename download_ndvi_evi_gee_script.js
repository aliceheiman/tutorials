// Point Coordinate
var roi = ee.Geometry.Point([lon, lat]);

// Polygon (drawn and named on the map)
//var roi = area_of_interest;

// Cloud + shadow mask using SCL + MSK_CLDPRB
// Source: Spatial Thoughts / S2_SR_HARMONIZED catalog
// Uses SCL (Scene Classification Layer)
function maskCloudAndShadows(image) {
  var cloudProb = image.select('MSK_CLDPRB');
  var snowProb  = image.select('MSK_SNWPRB');
  var scl       = image.select('SCL');

  var cloud  = cloudProb.lt(5);
  var snow   = snowProb.lt(5);
  var shadow = scl.eq(3);   // cloud shadow
  var cirrus = scl.eq(10);  // cirrus

  var mask = cloud.and(snow)
                  .and(shadow.neq(1))
                  .and(cirrus.neq(1));
  return image.updateMask(mask);
}

// Add NDVI and EVI bands
// S2 bands: B8 = NIR (10m), B4 = Red (10m)
// Source: GEE official tutorial (developers.google.com/earth-engine/tutorials/tutorial_api_06)
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('ndvi');
  
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      NIR:  image.select('B8').divide(10000),
      RED:  image.select('B4').divide(10000),
      BLUE: image.select('B2').divide(10000)
    }).rename('evi');

  return image.addBands([ndvi, evi]);
}

// Load S2 SR Harmonized collection
var collection = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(roi)
  .filterDate('2020-01-01', '2025-12-31') // Update to your desired window
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .map(maskCloudAndShadows)
  .map(addIndices);

// Export to CSV (for Python/R analysis) ─────────────────────────────────
// Source: Spatial Thoughts
var filteredCol = collection.select(['ndvi', 'evi']) // Derived bands
  .filter(ee.Filter.bounds(roi));

var timeSeries = ee.FeatureCollection(filteredCol.map(function(image) {
  var stats = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: roi,
    scale: 10,
    maxPixels: 1e10
  });

  // If the pixel is masked (cloudy), set ndvi to -9999 sentinel value
  var ndvi = ee.List([stats.get('ndvi'), -9999])
               .reduce(ee.Reducer.firstNonNull());
               
  var evi  = ee.List([stats.get('evi'),  -9999])
               .reduce(ee.Reducer.firstNonNull());

  return ee.Feature(null, {
    'date': ee.Date(image.get('system:time_start')).format('YYYY-MM-dd'),
    'ndvi': ndvi,
    'evi':  evi
  });
}));

// Export to drive
Export.table.toDrive({
  collection: timeSeries,
  description: 'S2_NDVI_EVI',
  folder: 'earthengine',
  fileNamePrefix: 'S2_NDVI_EVI',
  fileFormat: 'CSV',
  selectors: ['date', 'ndvi', 'evi']
});

// Chart NDVI time series at point ───────────────────────────────────────
// Source: ui.Chart.image.series
var chart = ui.Chart.image.series({
  imageCollection: collection.select('ndvi'),
  region: roi,
  reducer: ee.Reducer.mean(),   // or .first(), .median()
  scale: 10                      // S2 native resolution
}).setOptions({
  title: 'Sentinel-2 NDVI Time Series',
  vAxis: {title: 'NDVI', minValue: -0.1, maxValue: 1.0},
  hAxis: {title: 'Date', format: 'MMM-YYYY', gridlines: {count: 12}},
  lineWidth: 1,
  pointSize: 3,
  interpolateNulls: true   // connects across masked/cloudy gaps visually
});
print(chart);
