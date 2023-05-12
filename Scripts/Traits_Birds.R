# Script: Elton traits and PCA diet score of birds. Fixed Effect for BRMS
# Author: Iago Ferreiro Arias
# Date: 20th September 2022

###############################################################################
library(dplyr)
library(ggplot2)
##############################################################################

#Import AVONET  BirdLife data 
AVONET <- read.csv("Data/Traits/AVONET/AVONET_BirdLife_Data.csv", sep=";")
AVONET$ID_Avibase <- AVONET$Avibase.ID1

#Import AVONET eBird data 
ebird <- read.csv("Data/Traits/AVONET/AVONET_eBird_Data.csv", sep=";")
ebird$eBird_Species <- ebird$Species2
ebird$ID_Avibase <- ebird$Avibase.ID2
ebird<- ebird %>% dplyr::select(ID_Avibase, eBird_Species)

AVONET <- left_join(AVONET, ebird, by="ID_Avibase") #merge nomenclature
rm(ebird) #clean environment

#Import crosswalk between BirdLife dataset and BirdTree (Jetz phylo)
crosswalk <- read.csv("Data/Traits/AVONET/BirdLife_BirdTree_crosswalk.csv", sep=";")
crosswalk <- crosswalk %>% dplyr::select(BirdLife_Species, BirdTree_Species)

# merge with AVONET dataset
AVONET <- AVONET %>% rename(BirdLife_Species=Species1) 
AVONET <- left_join(AVONET, crosswalk, by="BirdLife_Species")
rm(crosswalk) #clean environment

# Relation between Hand Wing Index and Range size (proxy of avian flight efficiency and dispersal ability)
hist(log10(AVONET$Range.Size))
hist(sqrt(AVONET$Hand.Wing.Index))
model2 <-lm(log(Range.Size) ~ Hand.Wing.Index, data=AVONET)
sjPlot::plot_model(model2, type="pred")
rm(model2) #clean environment

# Relation between Kipps Distance and foraging stratum
hist(log(AVONET$Kipps.Distance))
ggplot(AVONET) + geom_boxplot(aes(x=Primary.Lifestyle, y=log(Kipps.Distance)), notch = T) +
  theme_classic()

#Reduce AVONET dataset to only predictors that we are interested
AVONET$Family <- AVONET$Family1
AVONET$Order <- AVONET$Order1
AVONET$Body_Mass <- AVONET$Mass                                                     
AVONET$Kipps_Dist <- AVONET$Kipps.Distance
AVONET$HWI <- AVONET$Hand.Wing.Index
AVONET$Trophic_Level <- AVONET$Trophic.Level
AVONET$Lifestyle <- AVONET$Primary.Lifestyle
AVONET$Tarsus_Length <- AVONET$Tarsus.Length

AVONET <- AVONET %>% dplyr::select(ID_Avibase, Order, Family, BirdLife_Species, 
                            eBird_Species, BirdTree_Species, Body_Mass, 
                            Tarsus_Length, HWI, Kipps_Dist, 
                            Lifestyle, Trophic_Level)

# Merge AVONET data with Elton traits (import "Elton" dataframe by running PCA_Birds.R)

setdiff(AVONET$BirdTree_Species, Elton$Species)
setdiff(AVONET$BirdLife_Species, Elton$Species)
setdiff(AVONET$eBird_Species, Elton$Species)

# merge by BirdTree species names (only 2 species do not coincide)
# check if Pampa curvipennis and pampa excellens are in BirdLife and eBird dataset
"Pampa curvipennis" %in% AVONET$BirdTree_Species 
"Pampa excellens" %in% AVONET$BirdTree_Species

#check if there is another name for those species in Elton Traits database
setdiff(Elton$Species, AVONET$BirdTree_Species)
# Campylopterus curvipennis = Pampa curvipennis
# Campylopterus excellens = Pampa excellens

# Change names in Elton dataset and merge datasets by BirdTree species names
Elton[which(Elton$Species=="Campylopterus curvipennis"), "Species"] <- "Pampa curvipennis"
Elton[which(Elton$Species=="Campylopterus excellens"), "Species"] <- "Pampa excellens" 

#check differences again
setdiff(AVONET$BirdTree_Species, Elton$Species) # solved
Elton <- Elton %>% rename(BirdTree_Species = Species)
TRAITS <- left_join(AVONET, Elton, by="BirdTree_Species")

rm(list = c("AVONET", "Elton")) #clean environment

# Merge TRAITS with evolutionary traits (phylo traits dataframe from Phylo_Birds.R)

setdiff(TRAITS$BirdTree_Species, phylo_traits$Species)

# Change names in Elton dataset and merge datasets by BirdTree species names
phylo_traits[which(phylo_traits$Species=="Campylopterus curvipennis"), "Species"] <- "Pampa curvipennis"
phylo_traits[which(phylo_traits$Species=="Campylopterus excellens"), "Species"] <- "Pampa excellens" 
setdiff(TRAITS$BirdTree_Species, phylo_traits$Species) #solved

phylo_traits <- phylo_traits %>% rename(BirdTree_Species = Species)
TRAITS <- left_join(TRAITS, phylo_traits, by="BirdTree_Species")

#check if NAs are coincident with species found in BirdLife nomenclature
setdiff(phylo_traits$BirdTree_Species,TRAITS$BirdTree_Species) %in% 
  TRAITS[which(is.na(TRAITS$Evo_Distinct)), "BirdLife_Species"]

#check if NAs are coincident with species found in eBird nomenclature
setdiff(phylo_traits$BirdTree_Species,TRAITS$BirdTree_Species) %in% 
  TRAITS[which(is.na(TRAITS$Evo_Distinct)), "eBird_Species"]

rm(phylo_traits)

#Import sheffers data: wildlife trade

Sheffers <- read.csv("Data/Traits/TRADE/Sheffers_aav5327_tables10_rev.csv",sep=";")
Sheffers <- subset(Sheffers, Class=="Bird")

Sheffers$Traded <- "Yes"

#check if differences are due to names
setdiff(Sheffers$Species, TRAITS$BirdLife_Species) %in% TRAITS$BirdTree_Species
setdiff(Sheffers$Species, TRAITS$BirdLife_Species) %in% TRAITS$eBird_Species

Sheffers <- Sheffers%>% dplyr::select(Species, Traded) %>% rename(BirdLife_Species = Species)

TRAITS <- left_join(TRAITS, Sheffers, by="BirdLife_Species")
Sheffers <- Sheffers %>% rename(BirdTree_Species= BirdLife_Species)
TRAITS <- left_join(TRAITS, Sheffers, by="BirdTree_Species")
Sheffers <- Sheffers %>% rename(eBird_Species= BirdTree_Species)
TRAITS <- left_join(TRAITS, Sheffers, by="eBird_Species")
TRAITS$Traded <- ifelse(is.na(TRAITS$Traded.x) & is.na(TRAITS$Traded.y) & is.na(TRAITS$Traded), "No", "Yes")

TRAITS <- TRAITS %>% dplyr::select(!Traded.x) %>% dplyr::select(!Traded.y)
table(TRAITS$Traded)

levels(as.factor(TRAITS$Trophic_Level))


TRAITS$Trophic_Level <- ifelse(TRAITS$Trophic_Level=="Carnivore", "Carnivore",
                        ifelse(TRAITS$Trophic_Level=="Omnivore", "Omnivore", 
                      ifelse(TRAITS$Trophic_Level=="Herbivore", "Herbivore", 
                                                    "Carnivore")))




rm(Sheffers)

writexl::write_xlsx(TRAITS, "Results/Traits_Birds_Projection.xlsx")


