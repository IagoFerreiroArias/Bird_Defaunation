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

small_df <- data.table::fread("Results/Output/Map_Hunting/small_df.csv")

nonhunted1 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_1.csv")
nonhunted2 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_2.csv")
nonhunted3 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_3.csv")
nonhunted4 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_4.csv")
nonhunted5 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_5.csv")
nonhunted6 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_6.csv")
nonhunted7 <- data.table::fread("Results/Output/Map_Hunting/nonhunted_small_7.csv")

# Split data into 1000 parts
num_splits <- 1000
split_size <- nrow(small_df)/num_splits

# Initialize an empty list to store results
result_list <- list()

# Process each split
for (i in 1:num_splits) {
  start_row <- (i - 1) * split_size + 1
  end_row <- i * split_size
  
  # Extract a chunk of the combined data
  combined_chunk <- small_df[start_row:end_row, ]
  print(paste0("Dataframe splitted: ", i))
  
  # Merge with small_df
  combined_chunk <- merge(combined_chunk, nonhunted1, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted2, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted3, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted4, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted5, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted6, by = "XY", all.x = TRUE)
  combined_chunk <- merge(combined_chunk, nonhunted7, by = "XY", all.x = TRUE)
  
  print(paste0("Dataframe merged: ", i))
  # Calculate mean for all columns excluding "XY"
  combined_chunk[, media := rowMeans(.SD, na.rm = TRUE), .SDcols = -"XY"]
  print(paste0("Mean calculated: ", i))
  # Calculate DI
  combined_chunk[, DI := 1 - media]
  print(paste0("DI calculated: ", i))
  #Split coordinates in two
  combined_chunk[, c("x", "y") := tstrsplit(XY, "_", fixed = TRUE)] 
  print(paste0("Coordinates splitted: ", i))
  # Extract XY and DI columns
  combined_chunk <- combined_chunk[, .(x, y, DI)]
  print(paste0("X, Y and DI extracted: ", i))
  # Append to the result list
  result_list[[i]] <- combined_chunk
}

# Bind the results into one dataset
combined_map <- rbindlist(result_list)

# Write the combined result to a CSV file
fwrite(combined_map, "Results/Output/Map_Defaunation/small_defaunation_map.csv")

#define original projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#rasterize defaunation index
tic()
combined_map <- raster::rasterFromXYZ(combined_map, crs = wgs84, res = 0.008333333) # Approximately 1 km at the Equator
toc()

#define projection 
mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tic()
projectRaster(combined_map, crs = mollweide, overwrite = TRUE,
              filename = "Results/Output/Map_Defaunation/small_defaunation_map.tif")
toc()
