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

