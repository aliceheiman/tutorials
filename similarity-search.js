// ════════════════════════════════════════════════════════════════════════════
// AlphaEarth Change Detection via Cosine Similarity (2017 vs 2024)
// ════════════════════════════════════════════════════════════════════════════

var CONFIG = {
  embeddingCollection: 'GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL',
  scale:               10,
  // Cosine similarity threshold: pixels BELOW this are flagged as "changed"
  // -1 = total opposite, 0 = unrelated, 1 = identical
  // Start around 0.5 and tune to taste
  changeThreshold:     0.6,
  // Minimum area (sq meters) for a changed region to be kept as a polygon
  minAreaM2:           1000,
  // Squareness filter: ratio of bbox area to polygon area
  // A perfect square = 1.0; natural shapes are typically < 0.6
  squarenessMin:       0.5,
  output_folder:       'MAPBIOMAS-EXPORT',
  exportPrefix:        'alphaearth_change',
};

// ── LOAD AOI ────────────────────────────────────────

var aoi = geometry

// ── SENTINEL-2 TRUE COLOUR (visual QA) ───────────────────────
// Cloud-masked median composite for each year over the AOI

function getS2Composite(year) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filter(ee.Filter.calendarRange(year, year, 'year'))
    // Dry season only — adjust months to suit your region
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .map(function(img) {
      // Basic cloud mask using the SCL band
      var scl = img.select('SCL');
      var cloudMask = scl.neq(3)   // cloud shadow
                        .and(scl.neq(8))   // medium cloud
                        .and(scl.neq(9))   // high cloud
                        .and(scl.neq(10)); // cirrus
      return img.updateMask(cloudMask).divide(10000); // scale to [0,1]
    })
    .median()
    .select(['B4', 'B3', 'B2'])  // RGB
    .clip(aoi);
}

var s2_2017 = getS2Composite(2017);
var s2_2024 = getS2Composite(2024);

var s2Vis = { min: 0.0, max: 0.25, gamma: 1.4 };

Map.addLayer(s2_2017, s2Vis, 'Sentinel-2 2017 (true colour)');
Map.addLayer(s2_2024, s2Vis, 'Sentinel-2 2024 (true colour)');

// ── LOAD EMBEDDINGS FOR 2017 AND 2024 ────────────────────────

var col = ee.ImageCollection(CONFIG.embeddingCollection);

function getYearEmbedding(year) {
  return col
    .filter(ee.Filter.date(
      ee.Date.fromYMD(year, 1, 1),
      ee.Date.fromYMD(year, 12, 31)
    ))
    .filter(ee.Filter.bounds(aoi))
    .mosaic()
    .clip(aoi);
}

var emb2017 = getYearEmbedding(2017);
var emb2024 = getYearEmbedding(2024);

print('Embedding bands:', emb2017.bandNames());

// ── COSINE SIMILARITY ─────────────────────────────────────────
// cos(θ) = (A · B) / (|A| * |B|)
// AlphaEarth embeddings are already unit-length, so |A| = |B| = 1
// meaning cos(θ) = A · B  (simple dot product)
// But we compute the norms anyway for robustness

var cosineSimilarity = emb2017.multiply(emb2024).reduce(ee.Reducer.sum());
print('Cosine similarity image:', cosineSimilarity);
print(ui.Chart.image.histogram(cosineSimilarity, aoi, 30))

// ── VISUALISE SIMILARITY MAP ──────────────────────────────────
// Diverging palette: blue = similar (no change), red = dissimilar (high change)

var simVis = {
  min: -1, max: 1,
  palette: ['#d73027','#f46d43','#fdae61','#fee08b',
            '#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850'],
};

Map.addLayer(cosineSimilarity, simVis, 'Cosine Similarity 2017→2024');

// ── EXPORT SIMILARITY RASTER ──────────────────────────────────

Export.image.toDrive({
  image:          cosineSimilarity,
  description:    CONFIG.exportPrefix + '_similarity_map',
  folder:         CONFIG.output_folder,
  fileNamePrefix: CONFIG.exportPrefix + '_similarity_map',
  region:         aoi,
  scale:          CONFIG.scale,
  crs:            'EPSG:4326',
  maxPixels:      1e10,
});

// ── THRESHOLD: FLAG HIGH-CHANGE PIXELS ───────────────────────
// Low cosine similarity = large change in land cover embedding

var changeMask = cosineSimilarity.lt(CONFIG.changeThreshold);

Map.addLayer(
  changeMask.selfMask(),
  {palette: ['orange']},
  'Changed pixels (below threshold)'
);

// ── VECTORISE CHANGED REGIONS ─────────────────────────────────

var changeVectors = changeMask
  .selfMask()
  .reduceToVectors({
    reducer:        ee.Reducer.countEvery(),
    geometry:       aoi,
    scale:          CONFIG.scale,
    maxPixels:      1e10,
    geometryType:   'polygon',
    eightConnected: true,   // connect diagonally adjacent pixels too
    labelProperty:  'change',
    tileScale:      4,
  });

// Filter out tiny fragments
changeVectors = changeVectors.filter(
  ee.Filter.gte('count', CONFIG.minAreaM2 / (CONFIG.scale * CONFIG.scale))
);

print('Changed regions (after min-area filter):', changeVectors.size());
Map.addLayer(changeVectors, {color: 'orange'}, 'Changed regions (vectorised)');

// ── SQUARENESS FILTER ─────────────────────────────────────────
// Squareness = polygon area / bounding-box area
// Ranges from 0 (very irregular) to 1.0 (perfect rectangle/square)
// Agricultural clearings and deforestation plots are typically > 0.7

var squareRegions = changeVectors.map(function(f) {
  var polyArea = f.geometry().area({maxError: 1});
  var bboxArea = f.geometry().bounds({maxError: 1}).area({maxError: 1});
  var squareness = polyArea.divide(bboxArea);
  return f.set('squareness', squareness).set('poly_area_m2', polyArea);
});

squareRegions = squareRegions.filter(
  ee.Filter.gte('squareness', CONFIG.squarenessMin)
);

print('Square-ish changed regions:', squareRegions.size());
Map.addLayer(squareRegions, {color: 'red'}, 'Square changed regions');

// Draw filled regions at low opacity, then overlay sharp outline
var filled = ee.Image().byte().paint(squareRegions, 1);
var outline = ee.Image().byte().paint(squareRegions, 1, 2); // last arg = stroke width

Map.addLayer(filled,   {palette: ['FF0000'], opacity: 0.15}, 'Square regions (fill)');
Map.addLayer(outline,  {palette: ['FF0000'], opacity: 0.9},  'Square regions (outline)');

// ── MAPBIOMAS LULC + SECONDARY VEG VISUALISATION (2017–2024) ─────────────

var lulcCollection = ee.Image(
  'projects/mapbiomas-public/assets/brazil/lulc/collection10_1/mapbiomas_brazil_collection10_1_coverage_v1'
);
var secVegCollection = ee.Image(
  'projects/mapbiomas-public/assets/brazil/lulc/collection10/mapbiomas_brazil_collection10_deforestation_secondary_vegetation_v2'
);

var lulcVisParams = {
  min: 0, max: 75,
  palette: [
    'ffffff','32a65e','32a65e','1f8d49','7dc975','04381d','026975',
    '000000','000000','7a6c00','ad975a','519799','d6bc74','d89f5c',
    'FFFFB2','edde8e','000000','000000','f5b3c8','C27BA0','db7093',
    'ffefc3','db4d4f','ffa07a','d4271e','db4d4f','0000FF','000000',
    '000000','ffaa5f','9c0027','091077','fc8114','2532e4','93dfe6',
    '9065d0','d082de','000000','000000','f5b3c8','c71585','f54ca9',
    'cca0d4','dbd26b','807a40','e04cfa','d68fe2','9932cc','e6ccff',
    '02d659','ad5100','000000','000000','000000','000000','000000',
    '000000','CC66FF','FF6666','006400','8d9e8b','f5d5d5','ff69b4',
    'ebf8b5','000000','000000','91ff36','7dc975','e97a7a','0fffe3',
    '000000','000000','000000','000000','000000','c12100',
  ]
};

var secVegVisParams = {
  min: 1, max: 7,
  palette: ['fffbc2','006400','4ea376','e31a1c','94fc03','ffa500','cccccc']
};

// MapBiomas band names follow the pattern "classification_YYYY"
var years = [2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024];

years.forEach(function(year) {
  var bandName = 'classification_' + year;

  var lulcLayer = lulcCollection
    .select(bandName)
    .clip(aoi);

  var secVegLayer = secVegCollection
    .select(bandName)
    .clip(aoi);

  // Add both layers but default to hidden — toggle manually in layer panel
  // Set the first and last years visible for quick before/after comparison
  var visible = (year === 2017 || year === 2024);

  //Map.addLayer(lulcLayer,   lulcVisParams,   'LULC '   + year, false);
  //Map.addLayer(secVegLayer, secVegVisParams, 'SecVeg ' + year, false);
});

Map.addLayer(lulcCollection, null, 'LULC', false);
Map.addLayer(secVegCollection, null, 'SecVeg', false);




// ── EXPORT VECTORS ────────────────────────────────────────────

// All changed regions
Export.table.toDrive({
  collection:     changeVectors,
  description:    CONFIG.exportPrefix + '_changed_regions',
  folder:         CONFIG.output_folder,
  fileNamePrefix: CONFIG.exportPrefix + '_changed_regions',
  fileFormat:     'GeoJSON',
});

// Square-ish regions only
Export.table.toDrive({
  collection:     squareRegions,
  description:    CONFIG.exportPrefix + '_square_regions',
  folder:         CONFIG.output_folder,
  fileNamePrefix: CONFIG.exportPrefix + '_square_regions',
  fileFormat:     'GeoJSON',
  selectors:      ['squareness', 'poly_area_m2', 'count'],
});
