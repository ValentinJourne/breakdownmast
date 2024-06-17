install.packages("ecmwfr")
library(ecmwfr)

# Set your API key
wf_set_key(user = "111279", key = "d884fd35-f644-426f-99ef-023daac5f92a", service = "cds")
#use can be your email from EOBS
 
# Define the request function for a specific location
# i will extract data from this daily up https://cds.climate.copernicus.eu/cdsapp#!/software/app-c3s-daily-era5-statistics?tab=overview

download_era5_data <- function(year, lat, lon) {
  request <- list(
    dataset_short_name = "reanalysis-era5-single-levels",
    product_type = "reanalysis",
    variable = "2m_temperature",
    year = as.character(year),
    month = sprintf("%02d", 1:12),
    day = sprintf("%02d", 1:31),
    time = sprintf("%02d:00", 0:23),
    format = "netcdf",
    target = paste0("era5_daily_temperature_", year, "_lat", lat, "_lon", lon, ".nc"),
    area = c(lat, lon, lat, lon)  # Define the bounding box for a specific point
  )
  
  wf_request(request, user = "111279")
}

# Loop through the years and download data for a specific location
start_year <- 1979  # Change this to the starting year of your interest
end_year <- 2023    # Change this to the ending year of your interest
latitude <- 52.52   # Replace with your specific latitude
longitude <- 13.41  # Replace with your specific longitude

plan(multicore, workers = parallel::detectCores())

years <- start_year:end_year

#download each layers in parallel processing 
future_lapply(years, download_era5_data, lat = latitude, lon = longitude)

#ok next option is with GEE 
#https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_DAILY#bands
#or this oone but issue because too big https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR
#temperature is in K and not degree
// Import the ECMWF ERA5 Daily Aggregated dataset
var dataset = ee.ImageCollection('ECMWF/ERA5/DAILY');

// Define the region of interest (replace with your specific coordinates)
var point = ee.Geometry.Point([13.41, 52.52]); // Longitude, Latitude

// Filter the dataset for the desired date range
var startDate = '1979-01-01';
var endDate = '2023-12-31';

var filteredDataset = dataset.filterDate(startDate, endDate)
.select('mean_2m_air_temperature');

// Function to extract temperature data for a point
var extractTemperature = function(image) {
  var date = ee.Date(image.get('system:time_start')).format('YYYY-MM-dd');
  var temp = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: point,
    scale: 1000,
    bestEffort: true,
    maxPixels: 1e9
  }).get('mean_2m_air_temperature');
  
  // Check if temp is null or undefined
  var tempExists = temp !== null && temp !== undefined;
  
  return ee.Feature(null, { 
    'date': date, 
    'temperature': ee.Algorithms.If(tempExists, temp, ee.Number.NaN)
  });
};

// Apply the extraction function to the dataset
var temperatureData = filteredDataset.map(extractTemperature);

// Convert the collection to a list
var temperatureList = temperatureData.toList(temperatureData.size());

// Print the first few entries
print(temperatureList);

// Export the data to Google Drive
Export.table.toDrive({
  collection: ee.FeatureCollection(temperatureList),
  description: 'TemperatureData',
  fileFormat: 'CSV'
});
