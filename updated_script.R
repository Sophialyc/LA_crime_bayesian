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

# read in csv for crime point csv
crime <- st_read('data/Crime_Data_from_2020_to_Present.csv') %>%
  st_as_sf(., coords = c("LON", "LAT"),
           crs= 4326) %>%
  st_transform(., 26945)
# check CRS 
st_crs(crime)
summary(crime)


# read in the zip code shapefile for LA county, the column of "zipcode" is for joining the population data
lapd <-st_read('data/clipped_zipcode/clipped.shp') %>%
  st_transform(., 26945)
# check the crs
st_crs(lapd)
qtm(lapd)


#plot the crime in the city, it may take a long time as there are 699280 points
tmap_mode("plot")
tm_shape(lapd) +
  tm_polygons(col = NA, alpha = 0.5) +
  tm_shape(crime) +
  tm_dots(col = "blue")

# read in population of LA of the US census data
population <- read.csv('data/2010_Census_Populations_by_Zip_Code.csv') 

# identify common columns
view(population) # Zip.code
view(lapd) # ZIPCODE

# clean up the names to easier merge
population <- clean_names(population) #zip_code
lapd <- clean_names(lapd) #zipcode
colnames(population)
names(population)[names(population) == "zip_code"] <- "zipcode"
# check the column names to verify the changes
colnames(population)

# join the population to the lapd shapefile
lapd_pop <- merge(lapd, population,
                  by.x="zipcode", 
                  by.y="zipcode",
                  no.dups = TRUE)

# lastly, read in the justice equity need index (JENI)
JENI <- read.csv('data/Justice_Equity_Need_Index_(zip_code).csv')
JENI <- clean_names(JENI)
names(JENI)[names(JENI) =="zip"] <- "zipcode"

# we only need the inequality driver index, therefore we will only extract the column with the zipcode out
JENI_driver <- subset(JENI, select = c(zipcode, driverspctl))

# merge with the la_pop shapefile
lapd_pop_JENI <- merge(lapd_pop, JENI_driver,
                  by.x="zipcode", 
                  by.y="zipcode",
                  no.dups = TRUE)

# now we have population, points of crime, inequality index and the shapefile of city of LA
# we have to get the column of point of crimes in each polygon 
# total number of crime in each zipcode
lapd_pop_JENI$crime <- lengths(st_intersects(lapd_pop_JENI, crime))
# check tge columns!
view(lapd_pop_JENI)

# delete unwanted columns, such as the url, zip
la_crime =subset(lapd_pop_JENI, select = -c(zip,tooltip,nla_url)) 

# calculate the expected number of cases
la_crime$exectednum_crime <- round(expected(population = la_crime$total_population, cases = la_crime$crime, n.strata =1), 0)

# reordering the columns
la_crime <- la_crime[, c(2,1,3,4,5,6,7,8,9,11,12,10)]
view(la_crime)

# change 0 in expected number of crime to 1 because the log function doen't work with 0
la_crime[la_crime$exectednum_crime == 0, "exectednum_crime"] <- 1
view(la_crime)

# coerce the dataframe into spatial object
sp.object <- as(la_crime, "Spatial")

# needs to be coerced into a matrix object
adjacencyMatrix <- shape2mat(sp.object)
# we extract the components for the ICAR model
extractComponents <- prep_icar_data(adjacencyMatrix)

# extract some components in the extractComponenets
n <- as.numeric(extractComponents$group_size)
nod1 <- extractComponents$node1
nod2 <- extractComponents$node2
n_edges <- as.numeric(extractComponents$n_edges)

# create dataset to be compiled in Stan
y <- la_crime$crime
x <- la_crime$driverspctl
e <- la_crime$exectednum_crime

# put all components into a list object
stan.spatial.dataset <- list(N=n, N_edges=n_edges, node1=nod1, node2=nod2, Y=y, X=x, E=e)

# use stan() to compile saved script
icar_poisson_fit = stan("icar_poisson_model.stan", data=stan.spatial.dataset, iter=2000, chains=4, verbose = FALSE)

# remove that annoying scientific notation
options(scipen = 999)
summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))$summary

summary(sp.object)
