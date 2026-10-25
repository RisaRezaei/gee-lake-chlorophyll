/**
 * Lake Chlorophyll-a Monitor (Sentinel-2 NDCI + MODIS-Aqua)
 * ---------------------------------------------------------
 * Tracks chlorophyll-a in a lake or reservoir with the Normalized Difference
 * Chlorophyll Index (NDCI, Mishra & Mishra 2012) from Sentinel-2, and compares
 * it with the MODIS-Aqua ocean-colour chlorophyll product where available.
 *
 * Outputs
 *   - Chart: monthly lake median and 90th percentile of Chl-a (Sentinel-2)
 *     with the MODIS-Aqua monthly mean on the same axis.
 *   - Map: mean Chl-a of a chosen season, plus the water mask.
 *   - CSV: month, chl_median, chl_p90, ndci_median, modis_chl, valid_km2.
 *   - Optional GeoTIFF of the seasonal Chl-a map.
 *
 * Chl-a from NDCI uses the published polynomial
 *   Chl-a (mg/m3) = 14.039 + 86.115 * NDCI + 194.325 * NDCI^2
 * which was fitted for turbid productive waters. Treat values as an index of
 * bloom intensity unless you calibrate the coefficients with field data.
 *
 * Author: Parisa Abbasalizadehrezaei · https://risarezaei.github.io
 * License: MIT
 */

// =============================== CONFIG ===================================
var CONFIG = {
  // Lake area (example: western basin of Lake Erie). Land is masked out.
  lake: ee.Geometry.Rectangle([-83.48, 41.40, -82.70, 42.02]),

  // Interactive charts stay fast for a few years; for long records (e.g. from
  // 2019) rely on the CSV export, which runs as a batch task.
  start: '2023-01-01',
  end: '2026-01-01',
  months: [5, 10],              // only these calendar months (inclusive); [1, 12] = all

  minCloudScore: 0.65,          // Cloud Score+ cs_cdf threshold (0-1)
  minWaterOccurrence: 80,       // JRC occurrence (%) for the permanent-water mask

  // Polynomial Chl-a = c0 + c1*NDCI + c2*NDCI^2 (Mishra & Mishra, 2012)
  coeffs: [14.039, 86.115, 194.325],

  seasonMap: {year: 2025, months: [7, 9]},  // months inclusive
  reduceScale: 100,             // metres for the lake statistics
  includeModis: true,

  exportTable: true,
  exportSeasonMap: false,
  driveFolder: 'earthengine',
  filePrefix: 'lake_chla'
};
// ==========================================================================

var lake = CONFIG.lake;
var start = ee.Date(CONFIG.start), end = ee.Date(CONFIG.end);

// Permanent water: avoids shoreline mixing and seasonal flooding.
var water = ee.Image('JRC/GSW1_4/GlobalSurfaceWater').select('occurrence')
  .gte(CONFIG.minWaterOccurrence).selfMask().clip(lake);

// ---- Sentinel-2 surface reflectance with Cloud Score+ -------------------
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
// Linking Cloud Score+ is done per time window: joining the whole period at
// once is much slower.
function s2(from, to) {
  return ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(lake).filterDate(from, to)
    .linkCollection(csPlus, ['cs_cdf'])
    .map(function (img) {
      var clear = img.select('cs_cdf').gte(CONFIG.minCloudScore);
      var notIce = img.select('SCL').neq(11);
      var ndci = img.normalizedDifference(['B5', 'B4']).rename('NDCI');
      var c = CONFIG.coeffs;
      var chl = ndci.expression('c0 + c1 * x + c2 * x * x',
        {x: ndci, c0: c[0], c1: c[1], c2: c[2]}).max(0).rename('CHL');
      return ndci.addBands(chl).updateMask(clear.and(notIce)).updateMask(water)
        .copyProperties(img, ['system:time_start']);
    });
}

// ---- MODIS-Aqua chlorophyll (NASA OBPG, ~4.6 km) ------------------------
var modis = ee.ImageCollection('NASA/OCEANDATA/MODIS-Aqua/L3SMI')
  .filterDate(start, end).select('chlor_a');

// ---- Monthly statistics --------------------------------------------------
// Month list built on the client, keeping only the chosen calendar months.
var monthStarts = [];
for (var d = new Date(CONFIG.start + 'T00:00:00Z'); d < new Date(CONFIG.end + 'T00:00:00Z');
     d.setUTCMonth(d.getUTCMonth() + 1)) {
  var mo = d.getUTCMonth() + 1;
  if (mo >= CONFIG.months[0] && mo <= CONFIG.months[1]) monthStarts.push(d.getTime());
}
var pixKm2 = CONFIG.reduceScale * CONFIG.reduceScale / 1e6;

var monthlyFc = ee.FeatureCollection(ee.List(monthStarts).map(function (d) {
  var m0 = ee.Date(d), m1 = m0.advance(1, 'month');
  var sub = s2(m0, m1);
  var img = ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.median(),
    ee.Image.constant([0, 0]).rename(['NDCI', 'CHL']).updateMask(0)));
  // One pass: median, 90th percentile and the number of clear water pixels.
  var stats = img.reduceRegion({
    reducer: ee.Reducer.median().combine(ee.Reducer.percentile([90]), null, true)
      .combine(ee.Reducer.count(), null, true),
    geometry: lake, scale: CONFIG.reduceScale, maxPixels: 1e13, tileScale: 4
  });
  var f = ee.Feature(null, {
    month: m0.format('YYYY-MM'), t: m0.millis(), n_images: sub.size(),
    chl_median: stats.get('CHL_median'), chl_p90: stats.get('CHL_p90'),
    ndci_median: stats.get('NDCI_median'),
    valid_km2: ee.Number(stats.get('CHL_count')).multiply(pixKm2)
  });
  if (CONFIG.includeModis) {
    var mod = modis.filterDate(m0, m1);
    var mv = ee.Algorithms.If(mod.size().gt(0),
      mod.mean().updateMask(water).reduceRegion({
        reducer: ee.Reducer.mean(), geometry: lake, scale: 4638.3, maxPixels: 1e10
      }).get('chlor_a'), null);
    f = f.set('modis_chl', mv);
  }
  return f;
}));

// Months with very little clear water give unreliable medians.
var good = monthlyFc.filter(ee.Filter.gt('valid_km2', 50));

print(ui.Chart.feature.byFeature(good, 't',
    CONFIG.includeModis ? ['chl_median', 'chl_p90', 'modis_chl'] : ['chl_median', 'chl_p90'])
  .setChartType('LineChart').setOptions({
    title: 'Monthly chlorophyll-a (mg/m³)',
    hAxis: {title: 'Month', format: 'yyyy'}, vAxis: {title: 'Chl-a (mg/m³)'},
    series: {0: {color: '2e7d32', lineWidth: 2}, 1: {color: 'a5d6a7', lineWidth: 1},
             2: {color: '1e88e5', lineWidth: 1.5, lineDashStyle: [4, 3]}},
    pointSize: 2, interpolateNulls: false
  }));
print('Monthly statistics', monthlyFc);

// ---- Seasonal map --------------------------------------------------------
var sm = CONFIG.seasonMap;
var seasonImg = s2(ee.Date.fromYMD(sm.year, sm.months[0], 1),
                   ee.Date.fromYMD(sm.year, sm.months[1], 1).advance(1, 'month'))
  .select('CHL').mean().clip(lake);

var CHL_VIS = {min: 0, max: 40, palette: ['08306b', '2171b5', '4eb3d3', '7fcdbb',
                                           'c7e9b4', 'fee08b', 'fc8d59', 'd73027']};
Map.centerObject(lake, 10);
Map.addLayer(water, {palette: ['0b3d91']}, 'Permanent water mask', false);
Map.addLayer(seasonImg, CHL_VIS, 'Chl-a mean ' + sm.year + ' months ' + sm.months.join('-'));

// ---- Exports -------------------------------------------------------------
if (CONFIG.exportTable) {
  Export.table.toDrive({
    collection: monthlyFc, description: CONFIG.filePrefix + '_monthly',
    folder: CONFIG.driveFolder, fileFormat: 'CSV',
    selectors: ['month', 'n_images', 'valid_km2', 'chl_median', 'chl_p90', 'ndci_median', 'modis_chl']
  });
}
if (CONFIG.exportSeasonMap) {
  Export.image.toDrive({
    image: seasonImg.toFloat(), description: CONFIG.filePrefix + '_season_' + sm.year,
    folder: CONFIG.driveFolder, region: lake, scale: 20, maxPixels: 1e13
  });
}
