# Script: Multivariate Environmental Similarity Surfaces 
# Author: Iago Ferreiro-Arias
# Date: 28th July, 2023

library(dplyr)
library(raster)
library(ggplot2)
library(data.table)
library(sp)
library(tictoc)
library(dismo)

#import hurdle model
load("Results/hurdle_model.RData") 

RR_data <- hurdle_model$data

#import rasters of tropical forest distribution
forest <- raster("Data/Spatial/Tropical_Forest/TropicalForest.tif")
forest2 <- shapefile ("Data/Spatial/Tropical_Forest/TropicalForest.shp")

#define projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

# Import raster layer for each predictor and extract values for each location.
# Subsequently, transform and rescale predictors for model projection
# scale and center using mean and sd values: raster values - mean(predictor hurdle) /sd(predictor hurdle)
df_rescale <- read.csv("Results/Rescaling_values.csv") #import mean and sd values for each predictor

#Distance to hunter's access points (euclidean distance to settlements, km).
DistHunt <- raster("Data/Spatial/Projections/DistHunt.tif")
DistHunt_log <- log10(DistHunt + 0.1) #log-transform while avoiding log10(0)
DistHunt_log <- (DistHunt_log- df_rescale[df_rescale$Predictor=="DistHunt", "Mean"])/df_rescale[df_rescale$Predictor=="DistHunt", "SD"]
DistHunt <- mask(DistHunt, forest2)
DistHunt_log <-projectRaster(from=DistHunt_log, to=forest, crs=wsg84, res=res(forest))
rm(DistHunt) #clean environment

#Travel time to major cities (>50k inhabitants) for 2015
TravDist <- raster("Data/Spatial/Projections/TravDist.tif")
TravDist_log <- log10(TravDist + 0.1) #log-transform while avoiding log10(0). #NAs produced since were set as -9999 in raster
TravDist_log <- (TravDist_log - df_rescale[df_rescale$Predictor=="TravDist", "Mean"])/df_rescale[df_rescale$Predictor=="TravDist", "SD"]
TravDist <- mask(TravDist, forest2)
TravDist_log <-projectRaster(from=TravDist_log, to=forest, crs=wsg84, res=res(forest))
rm(TravDist) #clean environment

#Prevalence of stunting among children < 5 years old. Up-to-date values from DHS, UNICEF etc. 
Stunting <- raster("Data/Spatial/Projections/Stunting.tif")
Stunting_log <- log10(Stunting + 0.1) #log-transform
Stunting_log <- (Stunting_log- df_rescale[df_rescale$Predictor=="Stunting", "Mean"])/df_rescale[df_rescale$Predictor=="Stunting", "SD"]
Stunting <- mask(Stunting, forest2)
Stunting_log <- projectRaster(from=Stunting_log, to=forest, crs=wsg84, res=res(forest))
rm(Stunting) #clean environment

#Human population density for 2020 (GPW V4.11)
PopDens <- raster("Data/Spatial/Projections/PopDens.tif")
# missing values from GPW comes from parks/protected areas, militry districts, airport zones, or non-reported in censuses.
# to avoid problems during predictors we set NA=0 as Benítez-López et al 2019.
PopDens[is.na(PopDens[])] <- 0  #replace NAs with 0
PopDens_log <- log10(PopDens + 0.1) #log-transform
PopDens_log <- (PopDens_log- df_rescale[df_rescale$Predictor=="PopDens", "Mean"])/df_rescale[df_rescale$Predictor=="PopDens", "SD"]
PopDens <- mask(PopDens, forest2)
PopDens_log <- projectRaster(from=PopDens_log, to=forest, crs=wsg84, res=res(forest))
rm(PopDens) #clean environment

#Net Primary productivity.
files <- list.files("Data/Spatial/NPP", pattern = "\\.TIFF$", full.names = TRUE)
NPP <- stack(files)
rm(files) #clean enviroment
NPP_mean <- mean(NPP, na.rm = TRUE) #estimate mean values for all years
writeRaster(NPP_mean, "Data/Spatial/Projections/NPP_mean.tif", overwrite=TRUE)

NPP_log <- log10(NPP_mean + 0.1) #log-transform while avoiding log10(0)
NPP_log <- (NPP_log - df_rescale[df_rescale$Predictor=="NPP", "Mean"])/df_rescale[df_rescale$Predictor=="NPP", "SD"]
NPP <- mask(NPP, forest2)
NPP_log <-projectRaster(from=NPP_log, to=forest, crs=wsg84, res=res(forest))
rm(NPP) #clean environment

#Reserves
Reserve <- raster("Data/Spatial/Projections/Reserves.tif")
Reserve[is.na(Reserve[])] <- 0 #outside reserves = NA
Reserve <- mask(Reserve, forest2)
Reserve <- projectRaster(from=Reserve, to=forest, crs=wsg84, res=res(forest))

#Country. Random effect. Code number for each country
CountryNum <- raster("Data/Spatial/Projections/Countries.tif")
CountryNum[is.na(CountryNum[])] <- 0 #outside reserves = NA
CountryNum <- mask(CountryNum, forest2)
CountryNum <- projectRaster(from=CountryNum, to=forest, crs=wsg84, res=res(forest))
#spatial_data[which(is.na(spatial_data$CountryNum)),"CountryNum"] <- 0

# prepare data for MESS: training and testing datasets
training <- RR_data %>% dplyr::select(DistHunt_log, TravDist_log, Stunting_log,
                                      PopDens_log, NPP_log)

testing <- raster::stack(DistHunt_log, TravDist_log, Stunting_log,
                         PopDens_log, NPP_log)

predictors <- stack(DistHunt_log, TravDist_log, Stunting_log,
                    PopDens_log, NPP_log, Reserve, CountryNum)

writeRaster(predictors, "Data/Spatial/Stack_Predictors_Scaled.tif", overwrite=TRUE)

tic()
MESS <- dismo::mess(x=testing,v=training, full=TRUE, filename = "Results/MESS.tif",
                    overwrite=TRUE)
toc()

#Clean environment
rm(list=c("DistHunt_log", "TravDist_log", "Stunting_log",
          "PopDens_log", "NPP_mean", "NPP_log", "Reserve", "CountryNum", 
          "training", "testing", "predictors"))