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

#Run first part of the code in CESGA Supercomputer

################# SMALL SPECIES #######################

rast <- raster::raster("Results/Output/Map_Defaunation/small_defaunation_map.tif")
resolution <- res(rast)
x <- resolution[1]/1000 #meters to km
y <- resolution[2]/1000 #meters to km
area_km2 <- x * y #pixel area in km2. each row in the datatable is a pixel
area_km2 <- as.numeric(area_km2)
print(paste0("Area of each pixel in km2: ",area_km2))
rm(list=c("x", "y")) #clean environment

small_map <- data.table::fread("Results/Output/Map_Defaunation/small_defaunation_map.csv")
small_map[, Realm := ifelse(x >= -20.00000 & x <= 60.00000, "Afrotropical",
                               ifelse(x > 60.00000 & x <= 165.00000, "Indomalayan", "Neotropical"))]
small_map[, Category := ifelse(DI >= 0 & DI < 0.1, 0.1,
                                ifelse(DI >= 0.1 & DI < 0.2, 0.2,
                                ifelse(DI >= 0.2 & DI < 0.3, 0.3,
                                ifelse(DI >= 0.3 & DI < 0.4, 0.4,
                                ifelse(DI >= 0.4 & DI < 0.5, 0.5,
                                ifelse(DI >= 0.5 & DI < 0.6, 0.6,
                                ifelse(DI >= 0.6 & DI < 0.7, 0.7,
                                ifelse(DI >= 0.7 & DI < 0.8, 0.8,
                                ifelse(DI >= 0.8 & DI < 0.9, 0.9,
                                ifelse(DI >= 0.9 & DI <= 1, 1, NA))))))))))]

# Conteo general por categoría (añadiendo Pantropical)
pantropical_small <- small_map[, .(Count = .N), by = .(Category)]
pantropical_small[, Realm := "Pantropical"]

# Conteo por categoría y Realm (3 niveles)
realm_small <- small_map[, .(Count = .N), by = .(Category, Realm)]

small_counts <- rbindlist(list(pantropical_small, realm_small),use.names=TRUE)

small_counts$Count <- as.numeric(small_counts$Count)
small_counts$Area_KM2 <- small_counts$Count * area_km2
small_counts[, Species := "Small"]


################# MEDIUM SPECIES #######################


rast <- raster::raster("Results/Output/Map_Defaunation/medium_defaunation_map.tif")
resolution <- res(rast)
x <- resolution[1]/1000 #meters to km
y <- resolution[2]/1000 #meters to km
area_km2 <- x * y #pixel area in km2. each row in the datatable is a pixel
area_km2 <- as.numeric(area_km2)
rm(list=c("x", "y")) #clean environment

medium_map <- data.table::fread("Results/Output/Map_Defaunation/medium_defaunation_map.csv")
medium_map[, Realm := ifelse(x >= -20.00000 & x <= 60.00000, "Afrotropical",
                            ifelse(x > 60.00000 & x <= 165.00000, "Indomalayan", "Neotropical"))]

medium_map[, Category := ifelse(DI >= 0 & DI < 0.1, 0.1,
                         ifelse(DI >= 0.1 & DI < 0.2, 0.2,
                         ifelse(DI >= 0.2 & DI < 0.3, 0.3,
                         ifelse(DI >= 0.3 & DI < 0.4, 0.4,
                         ifelse(DI >= 0.4 & DI < 0.5, 0.5,
                         ifelse(DI >= 0.5 & DI < 0.6, 0.6,
                         ifelse(DI >= 0.6 & DI < 0.7, 0.7,
                         ifelse(DI >= 0.7 & DI < 0.8, 0.8,
                         ifelse(DI >= 0.8 & DI < 0.9, 0.9,
                         ifelse(DI >= 0.9 & DI <= 1, 1, NA))))))))))]


# Conteo general por categoría (añadiendo Pantropical)
pantropical_medium <- medium_map[, .(Count = .N), by = .(Category)]
pantropical_medium[, Realm := "Pantropical"]

# Conteo por categoría y Realm (3 niveles)
realm_medium <- medium_map[, .(Count = .N), by = .(Category, Realm)]

medium_counts <- rbindlist(list(pantropical_medium, realm_medium), use.names=TRUE)
medium_counts$Count <- as.numeric(medium_counts$Count)
medium_counts$Area_KM2 <- medium_counts$Count * area_km2
medium_counts[, Species := "Medium"]

################# LARGE SPECIES #######################


rast <- raster::raster("Results/Output/Map_Defaunation/large_defaunation_map.tif")
resolution <- res(rast)
x <- resolution[1]/1000 #meters to km
y <- resolution[2]/1000 #meters to km
area_km2 <- x * y #pixel area in km2. each row in the datatable is a pixel
area_km2 <- as.numeric(area_km2)
rm(list=c("x", "y")) #clean environment

large_map <- data.table::fread("Results/Output/Map_Defaunation/large_defaunation_map.csv")
large_map[, Realm := ifelse(x >= -20.00000 & x <= 60.00000, "Afrotropical",
                            ifelse(x > 60.00000 & x <= 165.00000, "Indomalayan", "Neotropical"))]
large_map[, Category := ifelse(DI >= 0 & DI < 0.1, 0.1,
                        ifelse(DI >= 0.1 & DI < 0.2, 0.2,
                        ifelse(DI >= 0.2 & DI < 0.3, 0.3,
                        ifelse(DI >= 0.3 & DI < 0.4, 0.4,
                        ifelse(DI >= 0.4 & DI < 0.5, 0.5,
                        ifelse(DI >= 0.5 & DI < 0.6, 0.6,
                        ifelse(DI >= 0.6 & DI < 0.7, 0.7,
                        ifelse(DI >= 0.7 & DI < 0.8, 0.8,
                        ifelse(DI >= 0.8 & DI < 0.9, 0.9,
                        ifelse(DI >= 0.9 & DI <= 1, 1, NA))))))))))]


# Conteo general por categoría (añadiendo Pantropical)
pantropical_large <- large_map[, .(Count = .N), by = .(Category)]
pantropical_large[, Realm := "Pantropical"]

# Conteo por categoría y Realm (3 niveles)
realm_large <- large_map[, .(Count = .N), by = .(Category, Realm)]

large_counts <- rbindlist(list(pantropical_large, realm_large),use.names=TRUE)
large_counts$Count <- as.numeric(large_counts$Count)
large_counts$Area_KM2 <- large_counts$Count * area_km2
large_counts[, Species := "Large"]

################### ALL SPECIES #################

#afrotropical realm: longitude values from -20.00000 to 60.00000
#indomalayan realm: longitude values from 60.0000 to 165.00000

rast <- raster::raster("Results/Output/Map_Defaunation/combined_defaunation_map.tif")
resolution <- res(rast)
x <- resolution[1]/1000 #meters to km
y <- resolution[2]/1000 #meters to km
area_km2 <- x * y #pixel area in km2. each row in the datatable is a pixel
area_km2 <- as.numeric(area_km2)
rm(list=c("x", "y")) #clean environment

combined_map <- data.table::fread("Results/Output/Map_Defaunation/combined_defaunation_map.csv")
combined_map[, Realm := ifelse(x >= -20.00000 & x <= 60.00000, "Afrotropical",
                               ifelse(x > 60.00000 & x <= 165.00000, "Indomalayan", "Neotropical"))]
combined_map[, Category := ifelse(DI >= 0 & DI < 0.1, 0.1,
                           ifelse(DI >= 0.1 & DI < 0.2, 0.2,
                           ifelse(DI >= 0.2 & DI < 0.3, 0.3,
                           ifelse(DI >= 0.3 & DI < 0.4, 0.4,
                           ifelse(DI >= 0.4 & DI < 0.5, 0.5,
                           ifelse(DI >= 0.5 & DI < 0.6, 0.6,
                           ifelse(DI >= 0.6 & DI < 0.7, 0.7,
                           ifelse(DI >= 0.7 & DI < 0.8, 0.8,
                           ifelse(DI >= 0.8 & DI < 0.9, 0.9,
                           ifelse(DI >= 0.9 & DI <= 1, 1, NA))))))))))]


# Conteo general por categoría (añadiendo Pantropical)
pantropical_counts <- combined_map[, .(Count = .N), by = .(Category)]
pantropical_counts[, Realm := "Pantropical"]

# Conteo por categoría y Realm (3 niveles)
realm_counts <- combined_map[, .(Count = .N), by = .(Category, Realm)]

combined_counts <- rbindlist(list(pantropical_counts, realm_counts),use.names=TRUE)
combined_counts$Count <- as.numeric(combined_counts$Count)
combined_counts$Area_KM2 <- combined_counts$Count * area_km2
combined_counts[, Species := "All"]

print(colnames(combined_counts))
print(colnames(small_counts))
print(colnames(medium_counts))
print(colnames(large_counts))

counts <- rbindlist(list(combined_counts, small_counts, medium_counts, large_counts),use.names=TRUE)
fwrite(counts, "Results/Output/Map_Defaunation/Defaunation_Area_Data.csv")

########## BARPLOTS OF DEFAUNATION AREA ###################

library(ggplot2)


#Import data to create barplot
Area_Defaunation <- read.csv("Results/Output/Map_Defaunation/Defaunation_Area_Data.csv")
Area_Defaunation <- Area_Defaunation[-which(is.na(Area_Defaunation$Category)),]
Area_Defaunation$Area_KM2 <- Area_Defaunation$Area_KM2/1000000
Area_Defaunation$Category <- as.factor(Area_Defaunation$Category)


#### Pantropical ####
Pantropical <- Area_Defaunation[which(Area_Defaunation$Realm=="Pantropical"),]

g1 <- ggplot(data=Pantropical, aes(x=Category,y=Area_KM2, fill=Category)) +  
  geom_bar(stat = "identity")+ facet_grid(Species~.) + 
  theme_bw() + xlab("Defaunation index") + ylab("Area (million km2)") +
  theme(legend.position="None")+ labs(title="Indomalayan") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(legend.position="None") + labs(title="Pantropical") +
  scale_fill_manual(values=c("#6877ac","#919dba","#ffffbf","#f4d691",
                             "#f3cc85","#fdae61","#ee8a4e",
                             "#e6673b","#c53526","#ab1417"))


#### Neotropical ####
Neotropical <- Area_Defaunation[which(Area_Defaunation$Realm=="Neotropical"),]
  
g2<- ggplot(data=Neotropical, aes(x=Category,y=Area_KM2, fill=Category)) +  
    geom_bar(stat = "identity")+ facet_grid(Species~.) + 
    theme_bw() + xlab("Defaunation index") + ylab("Area (million km2)") +
  theme(legend.position="None")+ labs(title="Indomalayan") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="None")+ labs(title="Neotropical")+
  scale_fill_manual(values=c("#6877ac","#919dba","#ffffbf","#f4d691",
                              "#f3cc85","#fdae61","#ee8a4e",
                              "#e6673b","#c53526","#ab1417"))

#### Afrotropical ####
Afrotropical <- Area_Defaunation[which(Area_Defaunation$Realm=="Afrotropical"),]
  
g3<- ggplot(data=Afrotropical, aes(x=Category,y=Area_KM2, fill=Category)) +  
    geom_bar(stat = "identity")+ facet_grid(Species~.) + 
    theme_bw() + xlab("Defaunation index") + ylab("Area million km2") +
  theme(legend.position="None")+ labs(title="Indomalayan") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position="None")+ labs(title="Afrotropical") +
  scale_fill_manual(values=c("#6877ac","#919dba","#ffffbf","#f4d691",
                              "#f3cc85","#fdae61","#ee8a4e",
                              "#e6673b","#c53526","#ab1417"))
  
#### Indomalayan #####
Indomalayan <- Area_Defaunation[which(Area_Defaunation$Realm=="Indomalayan"),]
  
  g4<- ggplot(data=Indomalayan, aes(x=Category,y=Area_KM2, fill=Category)) +  
    geom_bar(stat = "identity")+ facet_grid(Species~.) + 
    theme_bw() + xlab("Defaunation index") + ylab("Area million km2") +
    theme(legend.position="None")+ labs(title="Indomalayan") +
    theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#6877ac","#919dba","#ffffbf","#f4d691",
                              "#f3cc85","#fdae61","#ee8a4e",
                              "#e6673b","#c53526","#ab1417"))
  
#### Merge plots ####
ggpubr::ggarrange(g1,g2,g3,g4, nrow=1, ncol=4, widths = c(0.8,0.8,0.8,0.8))

  
  
# Summary Area for each realm  
Area_Defaunation <- read.csv("Results/Output/Map_Defaunation/Defaunation_Area_Data.csv")
Area_Defaunation <- Area_Defaunation[-which(is.na(Area_Defaunation$Category)),]
Area_Defaunation<- Area_Defaunation %>% arrange(Realm, Category, Species)
Area_Defaunation$Category <- as.numeric(Area_Defaunation$Category)
Area_Defaunation$Defaunation <- ifelse(Area_Defaunation$Category ==0.1, "Intact",
                                ifelse(Area_Defaunation$Category >0.1 & Area_Defaunation$Category<0.3, "Low",
                                ifelse(Area_Defaunation$Category >=0.7, "High", "Moderate")))

Area_Defaunation$Area_KM2 <- Area_Defaunation$Area_KM2/1000000
df1 <- Area_Defaunation %>% group_by(Realm, Defaunation,Species) %>% summarise(Extent=sum(Area_KM2))

writexl::write_xlsx(df1,"Results/Output/Map_Defaunation/Defaunation_Area_Summary.xlsx")


####################### STATISTICS OF DEFAUNATION AREA ###########################

statistics <- data.frame("Realm"=c(rep("Pantropical",4),rep("Afrotropical",4),
                                   rep("Indomalayan",4),rep("Neotropical",4)),
                         "Species"=c(rep(c("All","Large", "Medium", "Small"),4)),
                         "Mean"=rep(NA,16),"Median"=rep(NA,16))
statistics$SD <- NA
statistics$IQR <- NA

mollweide  <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

all <- raster::raster("Results/Output/Map_Defaunation/combined_defaunation_map.tif")
large <- raster::raster("Results/Output/Map_Defaunation/large_defaunation_map.tif")
medium <-  raster::raster("Results/Output/Map_Defaunation/medium_defaunation_map.tif")
small <- raster::raster("Results/Output/Map_Defaunation/small_defaunation_map.tif")

statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="All"),"Mean"] <- mean(all[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="All"),"Median"] <- median(all[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="All"),"SD"] <- sd(all[], na.rm=TRUE)

statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Large"),"Mean"] <- mean(large[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Large"),"Median"] <- median(large[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Large"),"SD"] <- sd(large[], na.rm=TRUE)

statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Medium"),"Mean"] <- mean(medium[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Medium"),"Median"] <- median(medium[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Medium"),"SD"] <- sd(medium[], na.rm=TRUE)

statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Small"),"Mean"] <- mean(small[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Small"),"Median"] <- median(small[], na.rm=TRUE)
statistics[which(statistics$Realm=="Pantropical" & statistics$Species=="Small"),"SD"] <- sd(small[], na.rm=TRUE)

quantile(all[], na.rm=TRUE)

Africa <- sf::st_read("Data/Spatial/World_Continents/AfrotropicalRealm.shp")
Africa <- sf::st_transform(Africa, mollweide)
All_Africa <- raster::crop(all, Africa)
mean(All_Africa[], na.rm=TRUE)
gc()
median(All_Africa[], na.rm=TRUE)
gc()
sd(All_Africa[], na.rm=TRUE)
gc()
q<- quantile(All_Africa[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]


large_Africa <- raster::crop(large, Africa)
medium_Africa <- raster::crop(medium, Africa)
small_Africa <- raster::crop(small, Africa)

statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="All"),"Mean"] <- mean(All_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="All"),"Median"] <- median(All_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="All"),"SD"] <- sd(All_Africa[], na.rm=TRUE)

statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Large"),"Mean"] <- mean(large_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Large"),"Median"] <- median(large_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Large"),"SD"] <- sd(large_Africa[], na.rm=TRUE)


statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Medium"),"Mean"] <- mean(medium_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Medium"),"Median"] <- median(medium_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Medium"),"SD"] <- sd(medium_Africa[], na.rm=TRUE)

statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Small"),"Mean"] <- mean(small_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Small"),"Median"] <- median(small_Africa[], na.rm=TRUE)
statistics[which(statistics$Realm=="Afrotropical" & statistics$Species=="Small"),"SD"] <- sd(small_Africa[], na.rm=TRUE)


Asia <- sf::st_read("Data/Spatial/World_Continents/IndomalayanRealm.shp")
Asia <- sf::st_transform(Asia, mollweide)

All_Asia <- raster::crop(all, Asia)
All_Asia
All_Asia2 <- raster::mask(all,Asia)


mean(All_Asia[], na.rm=TRUE)
gc()
median(All_Asia[], na.rm=TRUE)
gc()
sd(All_Asia[], na.rm=TRUE)
gc()
q<- quantile(All_Asia[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]

gc()
All_Asia <- data.table::as.data.table(All_Asia, long = TRUE)
x <- terra::res(all)[1] / 1000
y <- terra::res(all)[2] / 1000
area_km2 <- x * y
sum(All_Asia$All_Asia >= 0.7) * area_km2
nrow(All_Asia) * area_km2

large_Asia <- raster::crop(large, Asia)
medium_Asia <- raster::crop(medium, Asia)
small_Asia <- raster::crop(small, Asia)

statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="All"),"Mean"] <- mean(All_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="All"),"Median"] <- median(All_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="All"),"SD"] <- sd(All_Asia[], na.rm=TRUE)


statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Large"),"Mean"] <- mean(large_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Large"),"Median"] <- median(large_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Large"),"SD"] <- sd(large_Asia[], na.rm=TRUE)


statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Medium"),"Mean"] <- mean(medium_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Medium"),"Median"] <- median(medium_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Medium"),"SD"] <- sd(medium_Asia[], na.rm=TRUE)

statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Small"),"Mean"] <- mean(small_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Small"),"Median"] <- median(small_Asia[], na.rm=TRUE)
statistics[which(statistics$Realm=="Indomalayan" & statistics$Species=="Small"),"SD"] <- sd(small_Asia[], na.rm=TRUE)

America <- sf::st_read("Data/Spatial/World_Continents/NeotropicalRealm.shp")
America <- sf::st_transform(America, mollweide)

All_America <- raster::crop(all, America)
mean(All_America[], na.rm=TRUE)
gc()
median(All_America[], na.rm=TRUE)
gc()
sd(All_America[], na.rm=TRUE)
gc()
q<- quantile(All_America[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]

large_America <- raster::crop(large, America)
medium_America <- raster::crop(medium, America)
small_America <- raster::crop(small, America)

statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="All"),"Mean"] <- mean(All_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="All"),"Median"] <- median(All_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="All"),"SD"] <- sd(All_America[], na.rm=TRUE)

statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Large"),"Mean"] <- mean(large_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Large"),"Median"] <- median(large_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Large"),"SD"] <- sd(large_America[], na.rm=TRUE)

statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Medium"),"Mean"] <- mean(medium_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Medium"),"Median"] <- median(medium_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Medium"),"SD"] <- sd(medium_America[], na.rm=TRUE)

statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Small"),"Mean"] <- mean(small_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Small"),"Median"] <- median(small_America[], na.rm=TRUE)
statistics[which(statistics$Realm=="Neotropical" & statistics$Species=="Small"),"SD"] <- sd(small_America[], na.rm=TRUE)
summary(small_America[])


#statistics of defaunation for large species
large <- raster::raster("Results/Output/Map_Defaunation/large_defaunation_map.tif")
mean(large[], na.rm=TRUE)
median(large[], na.rm=TRUE)
sd(large[], na.rm=TRUE)
q<- quantile(large[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]
rm(large)

#statistics of defaunation for medium species
medium <-  raster::raster("Results/Output/Map_Defaunation/medium_defaunation_map.tif")
mean(medium[], na.rm=TRUE)
median(medium[], na.rm=TRUE)
sd(medium[], na.rm=TRUE)
q<- quantile(medium[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]
rm(medium)
gc()

#statistics of defaunation for small species
small <- raster::raster("Results/Output/Map_Defaunation/small_defaunation_map.tif")
mean(small[], na.rm=TRUE)
median(small[], na.rm=TRUE)
sd(small[], na.rm=TRUE)
q<- quantile(small[], probs = seq(0, 1, 0.25), na.rm=TRUE)
q[4]-q[2]
