// ============================================================================
// Antarctic Minimum Snow Cover (MODIS MOD10A1 v061) - Clean Research Script
// Author: (edit) Your Name
// Last updated: 2026-02-23
//
// Purpose
// -------
// Compute the seasonal MINIMUM daily snow cover (%; NDSI_Snow_Cover) for Antarctica
// for each year, and export as GeoTIFFs to Google Drive.
//
// Dataset
// -------
// MODIS Terra Snow Cover Daily L3 Global 500 m, Version 6.1
// GEE ID: MODIS/061/MOD10A1
// Key bands used:
//   - NDSI_Snow_Cover                (0–100, percent; note: some fill values may exist)
//   - NDSI_Snow_Cover_Basic_QA       (QA categories; keep Best/Good/OK by default)
//
// Notes on QA
// -----------
// We use NDSI_Snow_Cover_Basic_QA to exclude low-quality retrievals.
// In practice, this band is commonly interpreted as categorical quality flags:
//   0 = Best, 1 = Good, 2 = OK, 3 = Other/poor (and/or fill depending on product)
// This script keeps QA <= 2 by default (Best/Good/OK).
// If you want a stricter mask, set KEEP_QA_MAX = 0 (Best only).
//
// Scientific window (editable)
// ----------------------------
// We compute MIN snow cover over an austral "melt season" window spanning year Y to Y+1.
// Default: Sep 23 (Y) to Mar 21 (Y+1).
// IMPORTANT: This window should be justified by your study design (phenology/seasonality).
//
// Exports
// -------
// Output images are exported at 500 m to Drive, in EPSG:3031 by default.
// You can change CRS/scale/region as needed.
// ============================================================================

// -------------------------------
// User parameters
// -------------------------------

// Region: Antarctica (south of 60°S). Adjust if needed.
var ANTARCTICA = ee.Geometry.Rectangle([-180, -90, 180, -60], null, false);

// Year range (inclusive)
var START_YEAR = 2002;
var END_YEAR   = 2023;

// Seasonal window definition (month/day) for start in year Y and end in year Y+1
var START_MONTH = 9;
var START_DAY   = 23;
var END_MONTH   = 3;
var END_DAY     = 21;

// QA threshold: keep QA values <= KEEP_QA_MAX
// 0 = Best only; 2 = Best/Good/OK
var KEEP_QA_MAX = 2;

// NoData fill value for exported rasters.
// Choose a value outside 0–100 (valid snow cover range). 255 is a common convention.
var NODATA_FILL = 255;

// Export settings
var EXPORT_FOLDER = 'GEE_exports';   // Drive folder name
var EXPORT_SCALE  = 500;            // meters
var EXPORT_MAXPIX = 1e13;           // allow large exports
var EXPORT_CRS    = 'EPSG:3031';    // Antarctic Polar Stereographic (commonly used)

// Optional: visualize one example year on the map
var PREVIEW_YEAR = 2015;
var VIS_PARAMS = {min: 0, max: 100, palette: ['000000', '0dffff', '0524ff', 'ffffff']};

// -------------------------------
// Helper functions
// -------------------------------

// Apply QA-based mask to a MOD10A1 image.
// Keeps pixels where NDSI_Snow_Cover_Basic_QA <= KEEP_QA_MAX.
function maskByBasicQA(img) {
  var qa = img.select('NDSI_Snow_Cover_Basic_QA');

  // Treat QA as categorical integer flags.
  // Keep Best/Good/OK by default (<= 2).
  var good = qa.lte(KEEP_QA_MAX);

  return img.updateMask(good);
}

// Compute seasonal minimum snow cover for year Y.
function seasonalMinSnowCover(year) {
  year = ee.Number(year);

  var startDate = ee.Date.fromYMD(year, START_MONTH, START_DAY);
  var endDate   = ee.Date.fromYMD(year.add(1), END_MONTH, END_DAY);

  var ic = ee.ImageCollection('MODIS/061/MOD10A1')
    .filterDate(startDate, endDate)
    .filterBounds(ANTARCTICA)
    .select(['NDSI_Snow_Cover', 'NDSI_Snow_Cover_Basic_QA'])
    .map(maskByBasicQA);

  // MIN over the season (daily product)
  var minSnow = ic.select('NDSI_Snow_Cover').min();

  // Fill masked pixels with NODATA_FILL to make missingness explicit in exports.
  minSnow = minSnow.unmask(NODATA_FILL).toUint8();

  // Add metadata
  minSnow = minSnow.set({
    'year': year,
    'window_start': startDate.format('YYYY-MM-dd'),
    'window_end': endDate.format('YYYY-MM-dd'),
    'keep_qa_max': KEEP_QA_MAX,
    'dataset': 'MODIS/061/MOD10A1'
  });

  return minSnow;
}

// Export a single year
function exportYear(year) {
  var img = seasonalMinSnowCover(year);

  var desc = 'MOD10A1_v061_Antarctica_MinSnowCover_Y' + year;

  Export.image.toDrive({
    image: img,
    description: desc,
    folder: EXPORT_FOLDER,
    fileNamePrefix: desc,
    region: ANTARCTICA,
    scale: EXPORT_SCALE,
    crs: EXPORT_CRS,
    maxPixels: EXPORT_MAXPIX
  });
}

// -------------------------------
// Main batch run
// -------------------------------

// Build a list of years and loop exports
var years = ee.List.sequence(START_YEAR, END_YEAR);
print('Exporting years:', years);
print('Season window:', START_MONTH + '/' + START_DAY + ' (Y) to ' + END_MONTH + '/' + END_DAY + ' (Y+1)');
print('QA threshold (keep <=):', KEEP_QA_MAX);
print('NoData fill value:', NODATA_FILL);

// Trigger exports (one task per year)
years.evaluate(function(listOfYears) {
  listOfYears.forEach(function(y) {
    print('Creating export task for year:', y);
    exportYear(y);
  });
});

// -------------------------------
// Optional preview
// -------------------------------
var preview = seasonalMinSnowCover(PREVIEW_YEAR);
Map.setCenter(0, -75, 3);
Map.addLayer(preview, VIS_PARAMS, 'Min snow cover (preview ' + PREVIEW_YEAR + ')');
Map.addLayer(ee.Image().paint(ANTARCTICA, 1, 2), {palette: ['FF0000']}, 'Antarctica ROI', false);
