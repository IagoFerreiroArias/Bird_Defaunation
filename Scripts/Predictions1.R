# Script: Spatial prediction of hunting impacts for tropical birds
# Author: Iago Ferreiro-Arias
# Date: 16th May 2023

#load libraries
library(raster)
library(sf)
library(tictoc)
library(dplyr)
library(data.table)
library(brms)
library(doParallel)
library(foreach)

# import data for model projections
spatial_data <- read.csv("Results/Spatial_data.csv") 
Mtraits <- read.csv("Results/Mtraits.csv") 

#check posible NAs in species names for phylogenetic tree
which(is.na(Mtraits$BirdTree_Species))
which(is.na(Mtraits$AOH_Species))
Mtraits <- Mtraits[-which(is.na(Mtraits$BirdTree_Species)),]

spatial_data$CountryNum <- as.factor (spatial_data$CountryNum)
spatial_data$Reserve <- as.factor (spatial_data$Reserve)

# Filter the dataframe for the target species
# Mtraits <- Mtraits[Mtraits$AOH_Species %in% target_species, ]
#rm(target_species)

#Try with 200 species
Mtraits <- Mtraits[1:200,]

#import hurdle model
load("Results/hurdle_model.RData") 

#import rasters of tropical forest distribution
forest <- raster("Data/Spatial/Tropical_Forest/TropicalForest.tif")

#define projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

### create cluster to run loop in parallel
print(parallel::detectCores()) #detect cores available for task (resources asked to CESGA)
n_cores <- parallel::detectCores() - 20 # avoid crash of memory in CESGA

cluster <- parallel::makeCluster(n_cores, typer="PSOCK") #define cluster so it can be used by %dopar%
# FORK backends do not work in clusters (?)
# PSOCK backends clonate environment for each core (not so efficient)=> Mantain environment clean during loop

doParallel::registerDoParallel(cl = cluster) #register it to be used by %dopar%
foreach::getDoParRegistered() #check if it is registered : TRUE?
print(foreach::getDoParWorkers()) #how many workers/cores are available?

tic("Entire loop") # elapsed time for entire loop

set.seed(12345) #reproducible analisis

foreach(i = Mtraits$AOH_Species) %dopar% {
  tictoc::tic(paste0("Running loop for",i))

  tictoc::tic(paste0("Prepare data for projections: ", i))
  
  proj_df <- data.frame(spatial_data, # dataframe with spatial predictors for projections
                        BirdTree_Species = rep(Mtraits[which(Mtraits$AOH_Species==i),"BirdTree_Species"], times = nrow(spatial_data)),
                        BodyMass_log = rep(Mtraits[which(Mtraits$AOH_Species==i),"BodyMass_log"], times = nrow(spatial_data)))
  
  print(paste0("Dataframe created for ",i))
  proj_df <- sp::SpatialPointsDataFrame(coords=cbind(proj_df$x,proj_df$y), 
                                        data=proj_df,proj4string=wgs84)
  print(paste0("Spatial Point Dataframe: ",i))
  
  #Charge AOH map for that species  (see Lumbierres et al 2022)
  AOH <- raster::raster(paste0("Data/Spatial/AOH_Birds/",i, ".tif"))
  
  # Calculate the factor to reduce the resolution: 1x1 km raster
  factor_x <- 0.008333333 / raster::res(AOH)[1] 
  factor_y <- 0.008333333 / raster::res(AOH)[2] 
  
  #set fun= "max" to return 1 if sp is present in at least 1 100x100m pixel (conservative)
  AOH_R <- raster::aggregate(AOH, fact=c(factor_x,factor_y), fun="max")
  print(paste0("AOH rasters created for ",i))
  
  #Extract values from AOH map to filter rows with AOH = NA. 
  proj_df$AOH <- raster::extract(AOH_R, proj_df, method="simple") #return values of the cell
  proj_df <- as.data.frame(proj_df) #return data as dataframe class
  
  print(paste0("AOH values extracted for ",i))
  
  if (all(is.na(proj_df$AOH))) {
    # If no point is within tropical forests, skip to the next species
    print(paste("Skipping species ", i, " as it is not present in tropical forests."))}
  else {
  print(paste("Predicting hunting impacts for ", i, " as it is present in tropical forests."))
  #Reduce size of dataframe for projections
  
  #proj_df <- proj_df[which(!is.na(proj_df$AOH)),] #filter rows with AOH = NA => Species not present
  proj_df <- na.omit(proj_df)
  tictoc::toc(log=TRUE)
  print(paste0("NA's in AOH removed for ",i))
  
  proj_df$DistHunt_log <- as.numeric(proj_df$DistHunt_log)
  proj_df$BodyMass_log <- as.numeric(proj_df$BodyMass_log)
  proj_df$TravDist_log <- as.numeric(proj_df$TravDist_log)
  proj_df$PopDens_log <- as.numeric(proj_df$PopDens_log)
  proj_df$Stunting_log <- as.numeric(proj_df$Stunting_log)
  proj_df$NPP_log <- as.numeric(proj_df$NPP_log)
  proj_df$Reserve <- as.factor(proj_df$Reserve)
  proj_df$CountryNum <- as.factor(proj_df$CountryNum)
  proj_df$BirdTree_Species <- as.factor(proj_df$BirdTree_Species)
  
  # perform predictions of hunting impacts for bird species
  tictoc::tic(paste0("Predicting hunting impacts for ", i))
  pred <- predict(hurdle_model, newdata= proj_df, allow_new_levels=TRUE, type= "response") #re_formula = NA)
  
  print(paste0("Starting preditction for ",i))
  #predicted values of RR.
  prediction <- data.frame(x = proj_df$x, y = proj_df$y, RR=round(pred[,1], digits=2))
  
  #truncate RR values in order to map only the declines in abundance RR<1
  pred_trunc<- prediction
  pred_trunc$RR <-ifelse(pred_trunc$RR > 1,1, pred_trunc$RR) # set(max(RR) = 1)
  tictoc::toc(log=TRUE)
  
  #rasterize predictions
  tictoc::tic(paste0("Mapping hunting impacts for ", i))
  raster_RR <- raster::rasterFromXYZ(prediction, crs = wgs84, res = 0.008333333) #1km raster
  #raster_RR <- raster::projectRaster(from=raster_RR, to=forest, crs=wgs84, res=res(forest))
  
  raster_RR_trunc <- raster::rasterFromXYZ(pred_trunc, crs = wgs84, res = 0.008333333)#1km raster
  #raster_RR_trunc <- raster::projectRaster(from=raster_RR_trunc, to=forest, crs=wgs84, res=res(forest))
  tictoc::toc(log=TRUE)
  
  #save predictions in 3 different folders.
  tictoc::tic(paste0("Saving",i, "raster")) #OUTPUT PER SPECIES
  raster::writeRaster(raster_RR, filename=paste0("Results/Output/Predictions/Hunted/pred_",i,".tif"),overwrite =TRUE) 
  raster::writeRaster(raster_RR_trunc,filename=paste0("Results/Output/Predictions_Truncated/Hunted/pred_trunc",i,".tif"),overwrite =TRUE) 
  tictoc::toc(log=TRUE)
  
  #clean environment
  rm(list=c("AOH", "factor_x","factor_y", "xy", "pred", "pred_trunc", 
            "proj_df", "pred", "prediction", "pred_trunc", "raster_RR",
            "raster_RR_trunc")) 
  
  warnings() #print warnings if any
  tictoc::toc(log = TRUE) #entire loop for a given bird species
}
}
tictoc::toc(log = TRUE) #entire loop