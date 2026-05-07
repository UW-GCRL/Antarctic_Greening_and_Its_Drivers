// Initialize Google Earth Engine region of interest (south of 60S)
var antarctica = ee.Geometry.Polygon(
  [[[-180, -60], [180, -60], [180, -90], [-180, -90]]], null, false
);

// Year of analysis (austral summer window: Sep 23 of `year` -> Mar 21 of `year+1`)
var year = 2004;

// Add NDVI and DOY (day-of-year) bands to each image
function calculateNdviWithDate(image) {
  var date = ee.Number.parse(image.date().format('D'));
  var red = image.select('Nadir_Reflectance_Band1').multiply(0.0001);
  var nir = image.select('Nadir_Reflectance_Band2').multiply(0.0001);
  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');
  return image.addBands(ndvi).addBands(ee.Image.constant(date).rename('Date').toInt());
}

// For each pixel, find the first and last DOY in the season where NDVI > 0.2
function findFirstAndLastValidDates(modis) {
  var firstDay = ee.Image().expression('9999').rename('FirstDay');
  var lastDay = ee.Image().expression('0').rename('LastDay');

  var dateImage = modis.iterate(function(image, previous) {
    image = ee.Image(image);
    var date = image.select('Date');
    var ndvi = image.select('NDVI').gt(0.2);

    var prev = ee.Dictionary(previous);
    var prevFirst = ee.Image(prev.get('FirstDay'));
    var prevLast = ee.Image(prev.get('LastDay'));

    var updatedFirst = prevFirst.where(ndvi.and(prevFirst.eq(9999)), date);
    var updatedLast = prevLast.where(ndvi, date);

    return ee.Dictionary({
      'FirstDay': updatedFirst,
      'LastDay': updatedLast
    });
  }, ee.Dictionary({'FirstDay': firstDay, 'LastDay': lastDay}));

  return {
    firstDay: ee.Image(ee.Dictionary(dateImage).get('FirstDay')),
    lastDay: ee.Image(ee.Dictionary(dateImage).get('LastDay'))
  };
}

// Build the MCD43A4 collection for the season and add NDVI / Date bands
var modis = ee.ImageCollection("MODIS/061/MCD43A4")
  .filterDate(ee.Date.fromYMD(year, 9, 23), ee.Date.fromYMD(year+1, 3, 21))
  .filterBounds(antarctica)
  .map(calculateNdviWithDate);

var results = findFirstAndLastValidDates(modis);
var firstValidDate = results.firstDay;
var lastValidDate = results.lastDay;

// Export helper: writes a single-band image to Drive at 500 m, MODIS sinusoidal
function exportImageToDrive(image, description, year) {
  Export.image.toDrive({
    image: image,
    description: description + '_' + year,
    scale: 500,
    region: antarctica,
    crs: 'SR-ORG:6974', // MODIS sinusoidal projection
    maxPixels: 1e9
  });
}


exportImageToDrive(firstValidDate, 'First_NDVI_above_02', year);
exportImageToDrive(lastValidDate, 'Last_NDVI_above_02', year);
