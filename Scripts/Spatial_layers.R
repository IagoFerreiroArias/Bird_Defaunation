#Script: Preparing spatial layers for predictions
#Author: Iago Ferreiro-Arias
#Date: 15th March, 2023

library(raster)
library(sp)
library(sf)
library(distanceto) #euclidean raster #install.packages("fasterize")
library(tictoc)
library(mailR)

#import mean and sd values for rescaling predictors
rescale <- read.csv("Data/Spatial/Rescaling_values.csv")

# Mask: distribution of tropical forests
forest <- shapefile("Data/Spatial/Tropical_Forest/TropicalForest.shp")

#### DISTANCE TO HUNTER ACCESS POINT #### 
tic()
#Import human settlement from WSF 2015
settle <- raster("Data/Spatial/HumanSettlement_WSF2015/WSF2015_v1_EPSG4326.vrt")
settle # check features and resolution
terra::distance(cbind(0,0), cbind(0,8.983046e-05), lonlat=T) # ~ 10 m resolution raster
settle <- mask(settle, forest)
toc()


# Euclidean distance to nearest human settlement
tic()
DistHunt <- distance(settle) #fill cells with NAs with euclidean distance
toc()
rm(settle) #clean environment

send.mail(from = "iago.ferreiro.arias@gmail.com",
          to = "iago.ferreiro.arias@gmail.com",
          subject = "Raster distance finished",
          body = "",
          smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "iago.ferreiro.arias@gmail.com", passwd = "OOJloErx", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#turn into ~1000m resolution raster
tic()
DistHunt <- raster::aggregate(DistHunt, fun=mean, fact=c(10,10)) # 10m x 10m (factor=100 times fewer cells)
toc()
DistHunt <- DistHunt/1000 #transform to km
DistHunt_log <- log10(DistHunt + 0.1) # log transform

# rescale raster with df values for predictions
DistHunt_log <- (DistHunt_log - rescale[which(rescale$Predictor=="DistHunt"), "Mean"])/rescale[which(rescale$Predictor=="DistHunt"), "SD"]
#save raster
writeRaster(DistHunt, filename = "Data/Spatial/HumanSettlement_WSF2015/DistHunt.tif", overwrite=T)
writeRaster(DistHunt_log, filename = "Data/Spatial/Projections/Scaled/DistHunt_log.tif", overwrite=T)            
rm(list=c("settle", "DistHunt", "DistHunt_log")) #clean environment
toc()

send.mail(from = "iago.ferreiro.arias@gmail.com",
          to = "iago.ferreiro.arias@gmail.com",
          subject = "Definitive raster distance saved",
          body = "",
          smtp = list(host.name = "smtp.gmail.com", port = 465, user.name = "iago.ferreiro.arias@gmail.com", passwd = "OOJloErx", ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#### TRAVEL DISTANCE TO NEAREST CITY #####
tic()
# import travel time to nearest cities
travel <- raster("Data/Spatial/TravelTime_2015/travel_time_to_cities_9.tif") # travel time to cities between 5k  and 10k inhabitants
travel
terra::distance(cbind(0,0), cbind(0,0.008333333), lonlat=T) #~ 1km resolution 
travel_log <- log10(travel + 0.1) #log transform    

# rescale raster with df values for predictions
TravDist_log <- (TravDist_log - rescale[which(rescale$Predictor=="TravDist"), "Mean"])/rescale[which(rescale$Predictor=="TravDist"), "SD"]

travel_t <- mask(TravDist_log, forest) # mask by forest distribution
writeRaster(travel_t, filename = "Data/Spatial/Projections/Scaled/TravDist_log.tif", overwrite=T) 
rm(list=c("travel", "travel_log", "TravelDist_log")) #clean environment 
toc()

#### PREVALENCE OF STUNTING ####
tic()
stunt <- shapefile("Data/Spatial/Stunting/Stunting_def.shp")
stunt <- raster::rasterize(stunt, travel_t,  field="stuntfinal", fun="mean")
stunt_log <- log10(stunt + 0.1) #log transform    
stunt_log <- (stunt_log - rescale[which(rescale$Predictor=="Stunting"), "Mean"])/rescale[which(rescale$Predictor=="Stunting"), "SD"]
stunt_log <- mask(sunt_log, forest)
writeRaster(stunt_log, filename== "Data/Spatial/Projections/Scaled/Stunting_log.tif", overwrite=T)
rm(list=c("stunt", "stunt_log")) #clean environment
toc()

#### COUNTRIES ####
tic()
countries <- shapefile("Data/Spatial/Countries/TropicalCountries.shp")
countries <- raster::rasterize(countries, travel_t,  field="CountryNum", fun="first") #rasterize
writeRaster(countries, filename = "Data/Spatial/Projections/Scaled/Countries.tif", overwrite=T) 
rm(countries) #clean environment
toc()


##### RESERVES #####
tic()
reserves <- shapefile("Data/Spatial/Reserves/Reserves_WDPA2022.shp")
reserves <- raster::rasterize(reserves, travel_t,  field="Reserve", background=0) #rasterize
reserves <- mask(reserves, forest)
writeRaster(reserves, filename = "Data/Spatial/Projections/Scaled/Reserves.tif", overwrite=T) 
toc()

#### FOOD BIOMASS ####
tic()
# upload density raster for different livestock species and multiply by average body mass
cattle <- raster("Data/Spatial/Food_Biomass/Cattle_km2.tif")
cattle <- cattle * 459.6 # cattle body mass

pig <- raster("Data/Spatial/Food_Biomass/Pig_km2.tif")
pig <- pig * 117.5  #pig body mass

chicken <- raster("Data/Spatial/Food_Biomass/Chicken_km2.tif")
chicken <- chicken * 2.112 # chicken body mass

goat <- raster("Data/Spatial/Food_Biomass/Goat_km2.tif")
goat <- goat * 54.9 #goat body mass

sheep<- raster("Data/Spatial/Food_Biomass/Sheep_km2.tif")
sheep <- sheep * 54.9 # sheep body mass

food_biomass <- cattle + pig + chicken + goat + sheep #aggregated livestock biomass raster

food_log <- log10(food_biomass + 0.1) #log transform            
# rescale raster with df values for predictions
FoodBiomass_log <- (food_log - rescale[which(rescale$Predictor=="FoodBiomass"), "Mean"])/rescale[which(rescale$Predictor=="FoodBiomass"), "SD"]

# mask by forest distribution and export raster
FoodBiomass_log <- mask(FoodBiomass_log,forest) 
writeRaster(FoodBiomass_log, filename = "Data/Spatial/Projections/Scaled/FoodBiomass_log.tif", overwrite=T) 
rm(list=c("cattle", "pig", "chicken", "goat", "sheep", "food_biomass", "food_log", "FoodBiomass_log")) # clean environment
toc()

#### HUMAN POPULATION DENSITY  ####
tic()
pop <- raster("Data/Spatial/GPW_v4_Pop_2020/GPW_v4_PopDens_2020.tif")          
pop_log <- log10(pop + 0.1)
PopDens_log <- (pop_log - rescale[which(rescale$Predictor=="PopDens"), "Mean"])/rescale[which(rescale$Predictor=="PopDens"), "SD"]

PopDens_log <- mask(PopDens_log, forest) #mask by tropical forest distribution
writeRaster(PopDens_log, filename = "Data/Spatial/Projections/Scaled/PopDens_log.tif", overwrite=T) #save raster
rm(list=c("pop", "pop_log", "PopDens_log", "forest")) #clean environment
toc()