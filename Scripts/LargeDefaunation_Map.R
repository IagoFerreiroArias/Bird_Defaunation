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

large_df <- data.table::fread("Results/Output/Map_Hunting/large_df.csv")
nonhunted_df <- data.table::fread("Results/Output/Map_Hunting/nonhunted_large_df.csv")

large <- merge(x=large_df, y=nonhunted_df, by = c("XY"), all.x=TRUE)
print("data merged")
print(ncol(large)-1)
large_map <-large[, media := rowMeans(.SD, na.rm = TRUE), .SDcols = -c("XY")]
#print(summary(large_map[,media]))
print("Mean calculated")

large_map <- large_map[, DI := 1 - media]
print("DI calculated")

large_map <- large_map[, c("x", "y") := tstrsplit(XY, "_", fixed = TRUE)] #split XY into two columns
#print("XY column explited into two columns")

large_map <- large_map[, .(x, y, DI)]
print("Columns x, y and DI selected")

fwrite(large_map, "Results/Output/Map_Defaunation/large_defaunation_map.csv")

#define original projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#rasterize defaunation index
tic()
large_map <- raster::rasterFromXYZ(large_map, crs = wgs84, res = 0.008333333) # Approximately 1 km at the Equator
toc()
#define projection and 

mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tic()
projectRaster(large_map, crs = mollweide, overwrite = TRUE,
              filename = "Results/Output/Map_Defaunation/large_defaunation_map.tif")
toc()