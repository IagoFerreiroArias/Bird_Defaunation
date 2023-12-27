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

# pet traded species
pet <- traits[-which(is.na(traits$Pet)),"AOH_Species"]
pet <- paste0(folder_path,"/pred_trunc_",pet,".csv")
pet_files <- csv_files[csv_files %in% pet] #select files of pet species
print(length(pet_files))
rm(pet) #clean environment

# Create an empty dataframe with the X and Y coordinates from spatial_data and round 5 decimals

#pet_df <- spatial_data %>% dplyr::select(x,y)
#pet_df$x <- round(pet_df$x ,5)
#pet_df$y <- round(pet_df$y ,5)
  
# Create a data.table from the merged_data for efficient merging
pet_df <- data.frame(XY=paste0(spatial_data$x,"_",spatial_data$y)) #create unique ID per grid cell
pet_df <- as.data.table(pet_df)

for(csv_file in pet_files) {
  csv <- data.table::fread(csv_file)  # Use fread for efficient reading
  csv <- csv[, .(x, y, RR)]
  #csv[, c("x", "y") := .(round(x, 5), round(y, 5))]  # Round coordinates
  colname <- gsub(".*/pred_trunc_|\\.csv", "", csv_file)
  
  # Set column names
  
  csv <- setnames(csv, new=c("x","y",colname))
  csv[, XY := paste(x, y, sep = "_")][, c("x", "y") := NULL]
  csv <- na.omit(csv)
  #print(all(csv$XY %in% pet_df$XY))
  
  # Debugging: Print data types and structure of csv
  print(paste0("Merging ", colname))
  
  #data.table::setkeyv(csv, c("x", "y"))  # Set key for merging
  #pet_df <- merge(x=pet_df, y=csv, by = c("x", "y"), all.x=TRUE)
  data.table::setkeyv(csv, c("XY"))  # Set key for merging
  pet_df <- merge(x=pet_df, y=csv, by = c("XY"), all.x=TRUE)
}
write.csv(pet_df,"Results/Output/Map_Hunting/pet_df.csv", row.names = FALSE)

################################################################################

#load libraries
library(raster)
library(sf)
library(tictoc)
library(dplyr)
library(doParallel)
library(foreach)
library(stringr)
library(data.table)

pet <- fread("Results/Output/Map_Hunting/pet_df.csv")
print(paste0("Number of pet species: ", ncol(pet)-2)) # number of species

pet_map <-pet[, media := rowMeans(.SD, na.rm = TRUE), .SDcols = -c("XY")]
print(summary(pet_map[,media]))
print("Mean calculated")

pet_map <- pet_map[, DI := 1 - media]
print("DI calculated")

pet_map <- pet_map[, c("x", "y") := tstrsplit(XY, "_", fixed = TRUE)] #split XY into two columns
#print("XY column explited into two columns")

pet_map <- pet_map[, .(x, y, DI)]
print("Columns x, y and DI selected")

fwrite(pet_map, "Results/Output/Map_Hunting/pet_map.csv")

#define original projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#rasterize defaunation index
tic()
pet_map <- raster::rasterFromXYZ(pet_map, crs = wgs84, res = 0.008333333) # Approximately 1 km at the Equator
toc()
#define projection and 

mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

tic()
projectRaster(pet_map, crs = mollweide, overwrite = TRUE,
              filename = "Results/Output/Map_Hunting/pet_map.tif")
toc()