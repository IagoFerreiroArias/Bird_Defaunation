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
library(ggplot2)

# import data for model projections
spatial_data <- read.csv("Results/Spatial_data.csv")
spatial_data$CountryNum <- as.factor(spatial_data$CountryNum)
spatial_data$Reserve <- as.factor(spatial_data$Reserve)

#Import bird species and traits
Mtraits <- read.csv("Results/Mtraits.csv") #hunted species filtered in Prep_Data.R
print(length(Mtraits$BirdTree_Species))

which(is.na(Mtraits$BirdTree_Species)) #check posible NAs in species names for phylogenetic tree
Mtraits <- Mtraits[-which(is.na(Mtraits$BirdTree_Species)),] #remove NAs
print(length(Mtraits$BirdTree_Species))

which(is.na(Mtraits$AOH_Species)) #check posible NAs in species names for rasters

# Supercomputer crash when we are trying to predict too large datasets in parallel.
# We will slipt the code of predictions for bird species with large range distributions
# (AOH files with large file size) and species with small range distributions
# To this end, we will use foreach package for small-sized files to run in paralel
# and for loops for large-sized files.

# Get a list of all AOH files in the directory
AOH_files <- list.files("Results/Output/Tropical_Hunted_AOH", pattern = "\\.tif$", full.names = TRUE)

# Extract species names from the AOH file paths (without ".tif" extension)
species_names <- tools::file_path_sans_ext(basename(AOH_files))

# filter dataframe for hunted species whose AOH is within tropical forest extension
Mtraits <- Mtraits[which(Mtraits$AOH_Species %in% species_names),]
print(length(Mtraits$BirdTree_Species))

# Get the file sizes corresponding to the species in Mtraits
file_sizes <- sapply(Mtraits$AOH_Species, function(species) {
  file <- file.path("Results/Output/Tropical_Hunted_AOH", paste0(species, ".tif"))
  file.info(file)$size / (1024^3) # File size in GB
})

# Add the "File_Size" column to Mtraits and sort it based on file size
Mtraits <- Mtraits %>% mutate(File_Size = file_sizes) %>% arrange(File_Size)
#write.csv(Mtraits, "Results/Birds_file_size.csv", row.names = FALSE)

Mtraits$BirdTree_Species <- as.factor(Mtraits$BirdTree_Species)

#import hurdle model
load("Results/hurdle_model.RData") 

#import rasters of tropical forest distribution
forest <- raster("Data/Spatial/Tropical_Forest/TropicalForest.tif")

#define projection for rasters
wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

### create cluster to run loop in parallel

tic("Entire loop") # elapsed time for entire loop

set.seed(12345) #reproducible analisis

#start projections
for(i in Mtraits$AOH_Species) {
  
  # Check if the output raster file already exists for the current species
  output_raster_file <- paste0("Results/Output/Predictions_Truncated/Hunted/pred_trunc_", i, ".tif")
  
  if (file.exists(output_raster_file)) {
    print(paste0("Raster for ", i, " already exists. Skipping prediction for this species."))}
  
  else{
    
    tictoc::tic(paste0("Prepare data for projections: ", i))
    
    proj_df <- data.frame(spatial_data, # dataframe with spatial predictors for projections from Prep_Data.R (predictors already scaled)
                          BirdTree_Species = rep(Mtraits[which(Mtraits$AOH_Species==i),"BirdTree_Species"], times = nrow(spatial_data)),
                          BodyMass_log = rep(Mtraits[which(Mtraits$AOH_Species==i),"BodyMass_log"], times = nrow(spatial_data)))
    
    print(paste0("Dataframe created for ",i))
    proj_df <- sp::SpatialPointsDataFrame(coords=cbind(proj_df$x,proj_df$y), 
                                          data=proj_df,proj4string=wgs84)
    print(paste0("Spatial Point Dataframe: ",i))
    
    #Charge rasampled and filtered AOH maps for hunted species
    AOH <- raster::raster(paste0("Results/Output/Tropical_Hunted_AOH/",i, ".tif"))
    
    #Extract values from AOH map to filter rows with AOH = NA. 
    proj_df$AOH <- raster::extract(AOH, proj_df, method="simple") #return values of the cell
    proj_df <- as.data.frame(proj_df) #return data as dataframe class
    
    print(paste0("AOH values extracted for ",i))
    
    if (all(is.na(proj_df$AOH))) {
      # If no point is within tropical forests, skip to the next species
      print(paste("Skipping species ", i, " as it is not present in tropical forests."))}
    else {
      print(paste("Predicting hunting impacts for ", i, " as it is present in tropical forests."))
      
      #Reduce size of dataframe for projections
      proj_df <- na.omit(proj_df) # omit sites not suitable for species and NA in predictors
      
      tictoc::toc(log=TRUE)
      print(paste0("NA's in AOH removed for ",i))
      
      if (nrow(proj_df) > 0 && all(!is.na(proj_df))) {
      #Make sure that predictors has same class as in brms model
      proj_df$DistHunt_log <- as.numeric(proj_df$DistHunt_log)
      proj_df$BodyMass_log <- as.numeric(proj_df$BodyMass_log)
      proj_df$TravDist_log <- as.numeric(proj_df$TravDist_log)
      proj_df$PopDens_log <- as.numeric(proj_df$PopDens_log)
      proj_df$Stunting_log <- as.numeric(proj_df$Stunting_log)
      proj_df$NPP_log <- as.numeric(proj_df$NPP_log)
      
      # Make sure factor levels are consistent with the original data
      proj_df$Reserve <- as.factor(proj_df$Reserve)
      proj_df$CountryNum <- as.factor(proj_df$CountryNum)
      proj_df$BirdTree_Species <- as.factor(proj_df$BirdTree_Species)
      
      # For large-sized files, predictions are too slow. Split dataset and run prediction in parallel
      
      # Divide dataframe in subdataset based on number of cores available
      n_cores <- strtoi(Sys.getenv("SLURM_CPUS_PER_TASK"))
      n_datasets <- 10
      subset_size <- ceiling(nrow(proj_df) / n_datasets)
      
      # Create a list with dataframe subsets
      subset_list <- split(proj_df, rep(1:n_datasets, each = subset_size, length.out = nrow(proj_df)))
      
      print(paste0("Split datasets created successfully for ",i))
      
      #Create empty list to save predictions
      predictions_list <- list()
      
      #create cluster
      cluster <- parallel::makeCluster(n_cores, typer="PSOCK") #define cluster so it can be used by %dopar%
      # FORK backends do not work in clusters (?)
      # PSOCK backends clonate environment for each core(?) (not so efficient)=> Mantain environment clean during loop
      
      doParallel::registerDoParallel(cl = cluster) #register it to be used by %dopar%
      foreach::getDoParRegistered() #check if it is registered : TRUE?
      print(foreach::getDoParWorkers()) #how many workers/cores are available?
      
      tic("Entire loop") # elapsed time for entire loop
      
      set.seed(12345) #reproducible analisis
      
      print(paste0("Starting prediction in parallel for ",i))
      
      output_dir <- "Results/Output/Predictions_Truncated/Hunted/Splitted/"
      
      predictions_list <- foreach(j = 1:length(subset_list)) %dopar% {
        splitted_df <- subset_list[[j]]
        
        # perform predictions of hunting impacts for bird species
        tictoc::tic(paste0("Predicting hunting impacts for ", i))
        pred <- predict(hurdle_model, newdata = splitted_df, allow.new.levels = TRUE, type = "response") #re_formula = NA != allow_new_levels
        
        # Store the results in a data frame
        result_df <- data.frame(x = splitted_df$x, y = splitted_df$y, RR = round(pred[, 1], digits = 2))
        
        # Generate the output file name
       
        file_number <- sprintf("%02d", j)
        output_file <- paste0(output_dir, "pred_trunc_", i, "_", file_number, ".csv")
        
        # Save the results as a CSV file with row.names = FALSE
        write.csv(result_df, file = output_file, row.names = FALSE)
        
        # Return any information you want to collect from this iteration
        # For example, you can return the output file name
        return(output_file)
      }
      
      # Stop the parallel processing
      stopCluster(cluster)
      print(paste0("Prediction in parallel performed for ",i))
      
      # list for storing dataframes
      csv_dataframes <- list()
      
      # CSV files with pattern of AOH_Species
      pattern <- paste0("pred_trunc_", i, "_.*\\.csv")
      matching_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
      
      # Upload files and save them into the list
      for (file in matching_files) {
        csv_dataframes[[length(csv_dataframes) + 1]] <- read.csv(file)
      }
      
      # combine dataframes into one
      prediction <- do.call(rbind, csv_dataframes)
      rm(list=c("csv_dataframes", "matching_files", "pattern"))
      print(str(prediction))
      print(paste0("Predictions merged for ",i))
      #prediction <- data.frame(x = proj_df$x, y = proj_df$y, RR=round(pred[,1], digits=2))
      
      #truncate RR values in order to map only the declines in abundance RR<1
      pred_trunc<- prediction
      print("pred_trunc created")
      pred_trunc$RR <-ifelse(pred_trunc$RR > 1,1, pred_trunc$RR) # set max(RR) = 1
      print("truncated RR vector created")
      tictoc::toc(log=TRUE)
      write.csv(pred_trunc,row.names = FALSE,
                file=paste0("Results/Output/Predictions_Truncated/Hunted/CSV/pred_trunc_",i,".csv")) 
      print(paste0("CSV saved for", i))
      
      #rasterize predictions
      tictoc::tic(paste0("Mapping hunting impacts for ", i))
      raster_RR <- raster::rasterFromXYZ(prediction, crs = wgs84, res = 0.008333333) #1km raster
      #raster_RR <- raster::projectRaster(from=raster_RR, to=forest, crs=wgs84, res=res(forest))
      print(paste0("Raster of predictions created for", i))
      
      raster_RR_trunc <- raster::rasterFromXYZ(pred_trunc, crs = wgs84, res = 0.008333333)#1km raster
      print(paste0("Raster of truncated predictions created for", i))
      #raster_RR_trunc <- raster::projectRaster(from=raster_RR_trunc, to=forest, crs=wgs84, res=res(forest))
      tictoc::toc(log=TRUE)
      
      #save predictions in 3 different folders.
      tictoc::tic(paste0("Saving",i, "raster")) #OUTPUT PER SPECIES
      raster::writeRaster(raster_RR, filename=paste0("Results/Output/Predictions/Hunted/pred_",i,".tif"),overwrite =TRUE) 
      raster::writeRaster(raster_RR_trunc,filename=paste0("Results/Output/Predictions_Truncated/Hunted/pred_trunc_",i,".tif"),overwrite =TRUE) 
      tictoc::toc(log=TRUE)
      
      #clean environment
      rm(list=c("AOH", "factor_x","factor_y", "xy", "pred", "pred_trunc", 
                "proj_df", "pred", "prediction", "pred_trunc", "raster_RR",
                "raster_RR_trunc")) 
      
      warnings() #print warnings if any
      tictoc::toc(log = TRUE) #entire task for a given bird species
      } 
      else {print(paste("Skipping predictions for species", i, "due to lack of suitable data."))}
    }
  }
}
tictoc::toc(log = TRUE) #entire parallel task

