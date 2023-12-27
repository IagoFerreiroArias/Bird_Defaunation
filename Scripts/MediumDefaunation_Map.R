# Script: Mapping hunting impacts in tropical birds
# Author: Iago Ferreiro-Arias
# Date: 10th August 2023

#load libraries
library(raster)
library(terra)
library(sf)
library(tictoc)
library(dplyr)
library(doParallel)
library(foreach)
library(stringr)
library(data.table)

medium_df <- data.table::fread("Results/Output/Map_Hunting/medium_df.csv")
nonhunted_df <- data.table::fread("Results/Output/Map_Hunting/nonhunted_medium_df.csv")

medium <- merge(x=medium_df, y=nonhunted_df, by = c("XY"), all.x=TRUE)
print("data merged")
print(ncol(medium)-1)
medium_map <-medium[, media := rowMeans(.SD, na.rm = TRUE), .SDcols = -c("XY")]
#print(summary(medium_map[,media]))
print("Mean calculated")

medium_map <- medium_map[, DI := 1 - media]
print("DI calculated")

medium_map <- medium_map[, c("x", "y") := tstrsplit(XY, "_", fixed = TRUE)] #split XY into two columns
#print("XY column explited into two columns")

medium_map <- medium_map[, .(x, y, DI)]
print("Columns x, y and DI selected")

fwrite(medium_map, "Results/Output/Map_Defaunation/medium_defaunation_map.csv")

#define original projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#rasterize defaunation index
tic()
medium_map <- raster::rasterFromXYZ(medium_map, crs = wgs84, res = 0.008333333) # Approximately 1 km at the Equator
toc()
#define projection and 

mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tic()
projectRaster(medium_map, crs = mollweide, overwrite = TRUE,
              filename = "Results/Output/Map_Defaunation/medium_defaunation_map.tif")
toc()