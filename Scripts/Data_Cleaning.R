<<<<<<< Updated upstream
# Script: Clean hunting database: Birds dataset
# Author: Iago Ferreiro Arias
# Date: 3th October 2022

###############################################################################
library(dplyr)
library(writexl) #export results into xlsx files
library(sp)
library(mapview)
##############################################################################


# Import hunting-induced defaunation data:
hunting_raw <- read.csv("Data/Hunting/Defaunation_db.csv", sep=";")

levels(as.factor(hunting_raw$Analisis))
hunting_raw <- subset(hunting_raw, Analisis=="yes") #remove rows not included in analisis

levels(as.factor(hunting_raw$Hloss_Fragm))
hunting_raw <- subset(hunting_raw, Hloss_Fragm=="None") #remove rows with potential confounding effects

levels(as.factor(hunting_raw$Group)) #filter bird species
hunting_birds <- subset(hunting_raw, Group=="Birds")
rm(hunting_raw) #clean environment

hunting_birds

#Create response variables
hunting_birds$RR <-hunting_birds$Hdens/hunting_birds$Udens # continuous models 
hunting_birds$bin <- ifelse(hunting_birds$Hdens == 0, 0, 1) # binomial models

#Create FoodBiomass predictor
hunting_birds$Cattle_biomass <- hunting_birds$CattleDens * 459.6
hunting_birds$Pig_biomass <- hunting_birds$PigDens * 117.5
hunting_birds$Chick_biomass<- hunting_birds$ChickDens * 2.112
hunting_birds$Goat_biomass <- hunting_birds$GoatDens * 54.9
hunting_birds$Sheep_biomass <- hunting_birds$SheepDens* 54.9

hunting_birds$FoodBiomass <- hunting_birds$Cattle_biomass + hunting_birds$Pig_biomass + 
  hunting_birds$Chick_biomass + hunting_birds$Goat_biomass + hunting_birds$Sheep_biomass 

#Select columns of interest

hunting_birds <- hunting_birds %>% dplyr::select(Study, Dataset=dataset2,
                                                Species, Nspp, RR, bin, 
                                                Hdens, Udens,
                                                Hunting=Type_of_hu, 
                                                Realm=Loc2,
                                                Country, Longitude, Latitude,
                                                Dist_Hunters = InterDist,
                                                TravDist, PopDens=PopDensUpd,
                                                Stunting= stunting2, FoodBiomass, 
                                                Reviewer, Reserve)


hunting_birds$Dist_Hunters <- as.numeric(hunting_birds$Dist_Hunters)/1000

table(hunting_birds$Nspp) #>200 ratio are grouped species (Eg. Mitu spp.)

hunting_birds <- subset(hunting_birds, Nspp==1)

#Import bird species taxonomy crosswalk
crosswalk <- read.csv("Data/Hunting/Crosswalk_Birds.csv", sep=";")
crosswalk <- crosswalk %>% dplyr::select(Species, BirdLife_Species, BirdTree_Species,
                                         IUCN_Species, IUCN_Status)

which(duplicated(crosswalk$Species))
hunting_birds <- left_join(hunting_birds, crosswalk, by="Species")
which(is.na(hunting_birds$IUCN_Status)) #check NAs

hunting_birds <- subset(hunting_birds, Udens!=0)

#Spatial overview
spatial <- SpatialPointsDataFrame(cbind(hunting_birds$Longitude, hunting_birds$Latitude),
                                  data=hunting_birds, 
                                  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

mapview(spatial, zcol="Reviewer", layer.name="Reviewer")
rm(spatial)


# Select specific traits (dataframe from Traits_Birds.R)

bird_traits <- TRAITS %>% dplyr::select(Order, Family, BirdTree_Species,
                                        Body_Mass, HWI, Kipps_Dist, Lifestyle,
                                        Trophic_Level, Activity, PCA_Dim1,
                                        PCA_Dim2, Evo_Distinctiveness,
                                        Tax_Distinctiveness, Traded)

bird_traits <- bird_traits[-which(duplicated(bird_traits$BirdTree_Species)),]#remove duplicated 
bird_traits[which(bird_traits$BirdTree_Species=="Stachyris chrysaea"),"BirdTree_Species"] <- "Cyanoderma chrysaeum"
#change name: Stachyris chrysaea == Cyanoderma chrysaeum
hunting_birds <- left_join(hunting_birds, bird_traits, by="BirdTree_Species")

which(is.na(hunting_birds$Traded)) #check join

hunting_birds$Realm <- ifelse(hunting_birds$Realm=="Africa", "Afrotropic",
                              ifelse(hunting_birds$Realm=="Asia", "Indomalayan", 
                                     "Neotropic"))

hunting_birds$Trophic_Level <- ifelse(hunting_birds$Trophic_Level=="Carnivore", "Carnivore",
                              ifelse(hunting_birds$Trophic_Level=="Omnivore", "Omnivore", 
                              ifelse(hunting_birds$Trophic_Level=="Herbivore", "Herbivore", 
                                     "Carnivore")))

hunting_birds <- hunting_birds %>% dplyr::select(Study, Dataset, Reviewer, Order, Family,
                                          Species, BirdLife_Species, BirdTree_Species,
                                          IUCN_Species, IUCN_Status,Hunting,RR, bin, 
                                          Hdens, Udens,Traded, Body_Mass, HWI, Kipps_Dist,
                                          Lifestyle, Trophic_Level, Activity,
                                          PCA_Dim1, PCA_Dim2, Evo_Distinctiveness,
                                          Tax_Distinctiveness, Realm,
                                          Country, Longitude, Latitude, Dist_Hunters,
                                          TravDist, PopDens, Stunting, FoodBiomass,
                                          Reserve)

length(unique(hunting_birds$Species))

hunting_birds$Reserve <- ifelse(hunting_birds$Reserve =="Yes", 1, 0)
levels(as.factor(hunting_birds$Reserve))

unique(hunting_birds$Country)

#numerical code for countries
countries <- read.csv("Data/Spatial/Countries.csv", sep=";")
unique(hunting_birds$Country) %in% countries$COUNTRY
countries <- countries %>% rename(Country= COUNTRY)
hunting_birds <- left_join(hunting_birds, countries, by="Country")
rm(countries)

#create dataset for modelling
RR_data <- hunting_birds %>% dplyr::filter(RR !="-Inf") %>% dplyr::filter(Hunting!="None")

#Spatial overview
spatial <- SpatialPointsDataFrame(cbind(RR_data$Longitude, RR_data$Latitude),
                                  data=RR_data, 
                                  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

mapview(spatial, zcol="Reviewer", layer.name="Reviewer")
rm(spatial)


rm(list=c("hunting_birds", "crosswalk", "bird_traits"))

write.csv(RR_data, "Data/Bird_RR_data.csv", row.names=F)



=======
# Script: Clean hunting database: Birds dataset
# Author: Iago Ferreiro Arias
# Date: 3th October 2022

###############################################################################
library(dplyr)
library(writexl) #export results into xlsx files
library(sp)
library(mapview)
##############################################################################


# Import hunting-induced defaunation data:
hunting_raw <- read.csv("Data/Hunting/Defaunation_db.csv", sep=";")

levels(as.factor(hunting_raw$Analisis))
hunting_raw <- subset(hunting_raw, Analisis=="yes") #remove rows not included in analisis

levels(as.factor(hunting_raw$Hloss_Fragm))
hunting_raw <- subset(hunting_raw, Hloss_Fragm=="None") #remove rows with potential confounding effects

levels(as.factor(hunting_raw$Group)) #filter bird species
hunting_birds <- subset(hunting_raw, Group=="Birds")
rm(hunting_raw) #clean environment

hunting_birds

#Create response variables
hunting_birds$RR <-hunting_birds$Hdens/hunting_birds$Udens # continuous models 
hunting_birds$bin <- ifelse(hunting_birds$Hdens == 0, 0, 1) # binomial models

#Create FoodBiomass predictor
hunting_birds$Cattle_biomass <- hunting_birds$CattleDens * 459.6
hunting_birds$Pig_biomass <- hunting_birds$PigDens * 117.5
hunting_birds$Chick_biomass<- hunting_birds$ChickDens * 2.112
hunting_birds$Goat_biomass <- hunting_birds$GoatDens * 54.9
hunting_birds$Sheep_biomass <- hunting_birds$SheepDens* 54.9

hunting_birds$FoodBiomass <- hunting_birds$Cattle_biomass + hunting_birds$Pig_biomass + 
  hunting_birds$Chick_biomass + hunting_birds$Goat_biomass + hunting_birds$Sheep_biomass 

#Select columns of interest

hunting_birds <- hunting_birds %>% dplyr::select(Study, Dataset=dataset2,
                                                Species, Nspp, RR, bin, 
                                                Hdens, Udens,
                                                Hunting=Type_of_hu, 
                                                Realm=Loc2,
                                                Country, Longitude, Latitude,
                                                Dist_Hunters = InterDist,
                                                TravDist, PopDens=PopDensUpd,
                                                Stunting= stunting2, FoodBiomass, 
                                                Reviewer, Reserve)


hunting_birds$Dist_Hunters <- as.numeric(hunting_birds$Dist_Hunters)/1000

table(hunting_birds$Nspp) #>200 ratio are grouped species (Eg. Mitu spp.)

hunting_birds <- subset(hunting_birds, Nspp==1)

#Import bird species taxonomy crosswalk
crosswalk <- read.csv("Data/Hunting/Crosswalk_Birds.csv", sep=";")
crosswalk <- crosswalk %>% dplyr::select(Species, BirdLife_Species, BirdTree_Species,
                                         IUCN_Species, IUCN_Status)

which(duplicated(crosswalk$Species))
hunting_birds <- left_join(hunting_birds, crosswalk, by="Species")
which(is.na(hunting_birds$IUCN_Status)) #check NAs

hunting_birds <- subset(hunting_birds, Udens!=0)

#Spatial overview
spatial <- SpatialPointsDataFrame(cbind(hunting_birds$Longitude, hunting_birds$Latitude),
                                  data=hunting_birds, 
                                  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

mapview(spatial, zcol="Reviewer", layer.name="Reviewer")
rm(spatial)


# Select specific traits (dataframe from Traits_Birds.R)

bird_traits <- TRAITS %>% dplyr::select(Order, Family, BirdTree_Species,
                                        Body_Mass, HWI, Kipps_Dist, Lifestyle,
                                        Trophic_Level, Activity, PCA_Dim1,
                                        PCA_Dim2, Evo_Distinctiveness,
                                        Tax_Distinctiveness, Traded)

bird_traits <- bird_traits[-which(duplicated(bird_traits$BirdTree_Species)),]#remove duplicated 
bird_traits[which(bird_traits$BirdTree_Species=="Stachyris chrysaea"),"BirdTree_Species"] <- "Cyanoderma chrysaeum"
#change name: Stachyris chrysaea == Cyanoderma chrysaeum
hunting_birds <- left_join(hunting_birds, bird_traits, by="BirdTree_Species")

which(is.na(hunting_birds$Traded)) #check join

hunting_birds$Realm <- ifelse(hunting_birds$Realm=="Africa", "Afrotropic",
                              ifelse(hunting_birds$Realm=="Asia", "Indomalayan", 
                                     "Neotropic"))

hunting_birds$Trophic_Level <- ifelse(hunting_birds$Trophic_Level=="Carnivore", "Carnivore",
                              ifelse(hunting_birds$Trophic_Level=="Omnivore", "Omnivore", 
                              ifelse(hunting_birds$Trophic_Level=="Herbivore", "Herbivore", 
                                     "Carnivore")))

hunting_birds <- hunting_birds %>% dplyr::select(Study, Dataset, Reviewer, Order, Family,
                                          Species, BirdLife_Species, BirdTree_Species,
                                          IUCN_Species, IUCN_Status,Hunting,RR, bin, 
                                          Hdens, Udens,Traded, Body_Mass, HWI, Kipps_Dist,
                                          Lifestyle, Trophic_Level, Activity,
                                          PCA_Dim1, PCA_Dim2, Evo_Distinctiveness,
                                          Tax_Distinctiveness, Realm,
                                          Country, Longitude, Latitude, Dist_Hunters,
                                          TravDist, PopDens, Stunting, FoodBiomass,
                                          Reserve)

length(unique(hunting_birds$Species))

hunting_birds$Reserve <- ifelse(hunting_birds$Reserve =="Yes", 1, 0)
levels(as.factor(hunting_birds$Reserve))

unique(hunting_birds$Country)

#numerical code for countries
countries <- read.csv("Data/Spatial/Countries.csv", sep=";")
unique(hunting_birds$Country) %in% countries$COUNTRY
countries <- countries %>% rename(Country= COUNTRY)
hunting_birds <- left_join(hunting_birds, countries, by="Country")
rm(countries)

#create dataset for modelling
RR_data <- hunting_birds %>% dplyr::filter(RR !="-Inf") %>% dplyr::filter(Hunting!="None")

#Spatial overview
spatial <- SpatialPointsDataFrame(cbind(RR_data$Longitude, RR_data$Latitude),
                                  data=RR_data, 
                                  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

mapview(spatial, zcol="Reviewer", layer.name="Reviewer")
rm(spatial)


rm(list=c("hunting_birds", "crosswalk", "bird_traits"))

write.csv(RR_data, "Data/Bird_RR_data.csv", row.names=F)



>>>>>>> Stashed changes
