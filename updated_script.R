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
icar_poisson_fit = stan("icar_poisson_model.stan", data=stan.spatial.dataset, iter=20000, chains=6, verbose = FALSE)

# remove that annoying scientific notation
options(scipen = 999)
# estimated result of the model
summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma"), probs=c(0.025, 0.975))$summary

# show first 6 rows only instead of the full 307
head(summary(icar_poisson_fit, pars=c("phi"), probs=c(0.025, 0.975))$summary)

print(icar_poisson_fit, pars=c("alpha", "beta", "sigma", "phi"), probs=c(0.025, 0.975))

# diagnostic check on the rHats - put everything into a data frame
diagnostic.checks <- as.data.frame(summary(icar_poisson_fit, pars=c("alpha", "beta", "sigma", "phi", "lp__"), probs=c(0.025, 0.5, 0.975))$summary)
# create binary variable
diagnostic.checks$valid <- ifelse(diagnostic.checks$Rhat < 1.1, 1, 0)
# tabulate it
table(diagnostic.checks$valid)

# show first 6 rows only instead of the full 307
head(summary(icar_poisson_fit, pars=c("mu"), probs=c(0.025, 0.975))$summary)

# extraction key posterior results for the generated quantities 
relativeRisk.results <- as.data.frame(summary(icar_poisson_fit, pars=c("mu"), probs=c(0.025, 0.975))$summary)
# now cleaning up this table up
# first, insert clean row numbers to new data frame
row.names(relativeRisk.results) <- 1:nrow(relativeRisk.results)
# second, rearrange the columns into order
relativeRisk.results <- relativeRisk.results[, c(1,4,5,7)]
# third, rename the columns appropriately
colnames(relativeRisk.results)[1] <- "rr"
colnames(relativeRisk.results)[2] <- "rrlower"
colnames(relativeRisk.results)[3] <- "rrupper"
colnames(relativeRisk.results)[4] <- "rHAT"

# view clean table 
head(relativeRisk.results)

# now, we proceed to generate our risk maps
# align the results to the areas in shapefile
la_crime$rr <- relativeRisk.results[, "rr"]
la_crime$rrlower <- relativeRisk.results[, "rrlower"]
la_crime$rrupper <- relativeRisk.results[, "rrupper"]

# create categories to define if an area has significant increase or decrease in risk, or nothing all 
la_crime$Significance <- NA
la_crime$Significance[la_crime$rrlower<1 & la_crime$rrupper>1] <- 0    # NOT SIGNIFICANT
la_crime$Significance[la_crime$rrlower==1 | la_crime$rrupper==1] <- 0  # NOT SIGNIFICANT
la_crime$Significance[la_crime$rrlower>1 & la_crime$rrupper>1] <- 1    # SIGNIFICANT INCREASE
la_crime$Significance[la_crime$rrlower<1 & la_crime$rrupper<1] <- -1   # SIGNIFICANT DECREASE

# For map design for the relative risk -- you want to understand or get a handle on what the distribution for risks look like
# this would inform you of how to create the labelling for the legends when make a map in tmap
summary(la_crime$rr)
hist(la_crime$rr)

# creating the labels
RiskCategorylist <- c(">0.0 to 0.25", "0.26 to 0.50", "0.51 to 0.75", "0.76 to 0.99", "1.00 & <1.01",
                      "1.01 to 1.10", "1.11 to 1.25", "1.26 to 1.50", "1.51 to 1.75", "1.76 to 2.00", "2.01 to 3.00")

# next, we are creating the discrete colour changes for my legends and want to use a divergent colour scheme
# scheme ranges from extreme dark blues to light blues to white to light reds to extreme dark reds
# you can pick your own colour choices by checking out this link [https://colorbrewer2.org]

RRPalette <- c("#65bafe","#98cffe","#cbe6fe","#dfeffe","white","#fed5d5","#fcbba1","#fc9272","#fb6a4a","#de2d26","#a50f15")

# categorising the risk values to match the labelling in RiskCategorylist object
la_crime$RelativeRiskCat <- NA
la_crime$RelativeRiskCat[la_crime$rr>= 0 & la_crime$rr <= 0.25] <- -4
la_crime$RelativeRiskCat[la_crime$rr> 0.25 & la_crime$rr <= 0.50] <- -3
la_crime$RelativeRiskCat[la_crime$rr> 0.50 & la_crime$rr <= 0.75] <- -2
la_crime$RelativeRiskCat[la_crime$rr> 0.75 & la_crime$rr < 1] <- -1
la_crime$RelativeRiskCat[la_crime$rr>= 1.00 & la_crime$rr < 1.01] <- 0
la_crime$RelativeRiskCat[la_crime$rr>= 1.01 & la_crime$rr <= 1.10] <- 1
la_crime$RelativeRiskCat[la_crime$rr> 1.10 & la_crime$rr <= 1.25] <- 2
la_crime$RelativeRiskCat[la_crime$rr> 1.25 & la_crime$rr <= 1.50] <- 3
la_crime$RelativeRiskCat[la_crime$rr> 1.50 & la_crime$rr <= 1.75] <- 4
la_crime$RelativeRiskCat[la_crime$rr> 1.75 & la_crime$rr <= 2.00] <- 5
la_crime$RelativeRiskCat[la_crime$rr> 2.00 & la_crime$rr <= 10] <- 6

# check to see if legend scheme is balanced - if a number is missing that categorisation is wrong!
table(la_crime$RelativeRiskCat)

# map of relative risk
rr_map <- tm_shape(la_crime) + 
  tm_fill("RelativeRiskCat", style = "cat", title = "Relavtive Risk", palette = RRPalette, labels = RiskCategorylist) +
  tm_shape(lapd) + tm_polygons(alpha = 0.05) + 
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.8, legend.text.size = 0.7) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))

# map of significance regions
sg_map <- tm_shape(la_crime) + 
  tm_fill("Significance", style = "cat", title = "Significance Categories", 
          palette = c("#33a6fe", "white", "#fe0000"), labels = c("Significantly low", "Not Significant", "Significantly high")) +
  tm_shape(lapd) + tm_polygons(col=NA, alpha = 0.10) +
  tm_layout(frame = FALSE, legend.outside = TRUE, legend.title.size = 0.9, legend.text.size = 0.6) +
  tm_compass(position = c("right", "top")) + tm_scale_bar(position = c("right", "bottom"))


# create side-by-side plot
tmap_arrange(rr_map, sg_map, ncol = 2, nrow = 1)
