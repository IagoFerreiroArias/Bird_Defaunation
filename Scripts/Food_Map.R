# Script: Mapping hunting impacts in tropical birds
# Author: Iago Ferreiro-Arias
# Date: 10th August 2023

#load libraries
library(raster)
library(sf)
library(tictoc)
library(dplyr)
library(doParallel)
library(foreach)
library(stringr)
library(data.table)

# Import the spatial_data and extract x, y columns
spatial_data <- read.csv("Results/Spatial_data.csv")

# Define the folder path
folder_path <- "Results/Output/Predictions_Truncated/Hunted/CSV"

# Get a list of CSV files in the folder
csv_files <- list.files(folder_path, pattern = "pred_trunc_.*\\.csv", full.names = TRUE)

# Import traits so select columns of interest
traits <- read.csv("Results/Traits_Birds_Projection.csv")

# food traded species
food <- traits[-which(is.na(traits$Food)),"AOH_Species"]
food <- paste0(folder_path,"/pred_trunc_",food,".csv")
food_files <- csv_files[csv_files %in% food] #select files of food species
print(paste0("Number of tropical species hunted for food:", length(food_files)))
rm(food) #clean environment

# Create an empty dataframe with the X and Y coordinates from spatial_data and round 5 decimals

#food_df <- spatial_data %>% dplyr::select(x,y)
#food_df$x <- round(food_df$x ,5)
#food_df$y <- round(food_df$y ,5)

# Create a data.table from the merged_data for efficient merging
food_df <- data.frame(XY=paste0(spatial_data$x,"_",spatial_data$y)) #create unique ID per grid cell
food_df <- as.data.table(food_df)

for(csv_file in food_files) {
  csv <- data.table::fread(csv_file)  # Use fread for efficient reading
  csv <- csv[, .(x, y, RR)]
  #csv[, c("x", "y") := .(round(x, 5), round(y, 5))]  # Round coordinates
  colname <- gsub(".*/pred_trunc_|\\.csv", "", csv_file)
  
  # Set column names
  
  csv <- setnames(csv, new=c("x","y",colname))
  csv[, XY := paste(x, y, sep = "_")][, c("x", "y") := NULL]
  csv <- na.omit(csv)
  #print(all(csv$XY %in% food_df$XY))
  
  # Debugging: Print data types and structure of csv
  print(paste0("Merging ", colname))
  
  #data.table::setkeyv(csv, c("x", "y"))  # Set key for merging
  #food_df <- merge(x=food_df, y=csv, by = c("x", "y"), all.x=TRUE)
  data.table::setkeyv(csv, c("XY"))  # Set key for merging
  food_df <- merge(x=food_df, y=csv, by = c("XY"), all.x=TRUE)
}

write.csv(food_df,"Results/Output/Map_Hunting/food_df.csv", row.names = FALSE)

food <- fread("Results/Output/Map_Hunting/food_df.csv")
food_map <-food[, media := rowMeans(.SD, na.rm = TRUE), .SDcols = -c("XY")]
print(summary(food_map[,media]))
print("Mean calculated")

food_map <- food_map[, DI := 1 - media]
print("DI calculated")

food_map <- food_map[, c("x", "y") := tstrsplit(XY, "_", fixed = TRUE)] #split XY into two columns
#print("XY column explited into two columns")

food_map <- food_map[, .(x, y, DI)]
print("Columns x, y and DI selected")

fwrite(food_map, "Results/Output/Map_Hunting/food_map.csv")

#define original projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#rasterize defaunation index
tic()
food_map <- raster::rasterFromXYZ(food_map, crs = wgs84, res = 0.008333333) # Approximately 1 km at the Equator
toc()
#define projection and 

mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tic()
projectRaster(food_map, crs = mollweide, overwrite = TRUE,
              filename = "Results/Output/Map_Hunting/food_map.tif")
toc()