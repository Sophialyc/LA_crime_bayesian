library(sf)
library(tmap)
library(tidyverse)
library(rgdal)
library(dplyr)
library(janitor)
library(readr)
library("spdep")
library("rstan")
library("geostan")
library("SpatialEpi")
library("tidybayes")
library(rgeos)


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)




# read in the shapefile for city of  LA
lapd <- st_read('data/LAPD/LAPD_Reporting_Districts.shp') %>%
  st_transform(., 26945)
# check CRS
st_crs(lapd)
# plot to see the shapefile
qtm(lapd)

# read in csv for crime point csv
crime <- st_read('data/Crime_Data_from_2020_to_Present.csv') %>%
  st_as_sf(., coords = c("LON", "LAT"),
           crs= 4326) %>%
  st_transform(., 26945)
# check CRS 
st_crs(crime)
summary(crime)


#plot the crime in the city, it may take a long time as there are 699280 points
tmap_mode("plot")
tm_shape(lapd) +
  tm_polygons(col = NA, alpha = 0.5) +
  tm_shape(crime) +
  tm_dots(col = "blue")

#dissolve the shapefile for uses later in cropping the zip_code boundary of LA
library(raster)
la_outline <- st_union(lapd, by_feature = FALSE)
ggplot(la_outline) + geom_sf(fill='grey')

outline <- st_boundary(la_outline)
qtm(outline)

# read in the zip code shapefile for LA county
la_zipcode <-st_read('data/Zip_Codes_(LA_County)/Zip_Codes_(LA_County).shp') %>%
  st_transform(., 26945)
# check the crs
st_crs(la_zipcode)
qtm(la_zipcode)

# since the shapefile is in LA county, which includes all other cities in LA country, hence, crop the area of city og LA using the outline created above
cropped <- st_intersection(la_zipcode, outline)
class(cropped)
qtm(cropped)
