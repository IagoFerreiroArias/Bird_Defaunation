library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(data.table)

# Cargar el archivo shapefile de países
countries <- st_read("Data/Spatial/Countries/TropicalCountries.shp")
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
countries <- st_transform(countries, mollweide)

# Cargar el raster

food_map <- rast("Results/Output/Map_Hunting/food_map.tif")

# Calcular las dimensiones de los píxeles en km
x <- res(food_map)[1] / 1000
y <- res(food_map)[2] / 1000

# Calcular el área de cada píxel en km^2
area_km2 <- x * y

#China <- crop(food_map, countries[which(countries$COUNTRY=="China"),])
#China <- data.table::as.data.table(China)
#sum(China$food_map>0.5, na.rm=TRUE) * area_km2

#Indonesia <- crop(food_map, countries[which(countries$COUNTRY=="Indonesia"),])
#Indonesia <- data.table::as.data.table(Indonesia)
#sum(Indonesia$food_map>0.5, na.rm=TRUE) * area_km2
#West_Africa <- countries %>% filter(COUNTRY=="Togo" |COUNTRY=="Nigeria" | 
#                                    COUNTRY=="Benin" | COUNTRY=="Ghana" |
#                                      COUNTRY=="Ivory Coast" | COUNTRY=="Guinea"|
#                                      COUNTRY=="Sierra Leone"|COUNTRY=="Guinea-Bissau")

#Africa<- crop(food_map,West_Africa)
#plot(Africa)
#Africa <- data.table::as.data.table(Africa)
#sum(Africa$food_map>0.5, na.rm=TRUE) * area_km2

# Iterar sobre cada país
result_list <- list()

for (i in 1:nrow(countries)) {
  # Obtener el polígono del país
  country_polygon <- st_geometry(countries[i, ])
  
  # Convertir el polígono a objeto terra
  terra_poly <- vect(st_as_sf(countries[i, ]))
  
  # Aplicar la máscara al raster usando el polígono
  masked_raster <- raster::mask(food_map, terra_poly)
  
  # Convertir el raster a data.table
  masked_dt <- as.data.table(masked_raster, long = TRUE)
  
  # Calcular el número total de píxeles forestales para el país
  forest_extent <- sum(!is.na(masked_dt$food_map), na.rm = TRUE)
  
  # Calcular el número de píxeles con DI > 0.3
  pixels_DI_gt_05 <- sum(masked_dt$food_map >= 0.3, na.rm = TRUE)
  
  # Calcular la columna Forest_Extent y Forest_km2
  country_df <- data.frame(
    Country = countries$COUNTRY[i],
    Forest_Extent = forest_extent,
    Forest_km2 = forest_extent * area_km2
  )
  
  # Calcular la variable DI
  DI <- pixels_DI_gt_05
  
  # Calcular DI_km2
  DI_km2 <- DI * area_km2
  
  # Calcular Perc_ForestDI05
  Perc_ForestDI05 <- (DI_km2 / country_df$Forest_km2) * 100
  
  # Agregar las nuevas columnas al data.frame
  country_df <- cbind(country_df, DI = DI, DI_km2 = DI_km2, Perc_ForestDI05 = Perc_ForestDI05)
  
  # Agregar la lista al resultado
  result_list[[i]] <- country_df
  gc()
}

# Combinar las listas en un solo data frame
result_df <- do.call(rbind, result_list)

# Ordenar el dataframe por DI_km2 en orden descendente
result_df <- result_df %>% arrange(desc(Perc_ForestDI05))
fwrite(result_df, "Results/Output/Map_Hunting/Food_Defaunation_Area.csv")

result_df <-fread("Results/Output/Map_Hunting/Food_Defaunation_Area.csv")
result_df <- result_df[-which(result_df$Country=="United States"),]
# Tomar los 35 países con mayores valores de DI_km2
top_food <- head(result_df, 35)
top_food$Country <- factor(top_food$Country, levels = top_food$Country[order(top_food$Perc_ForestDI05)])

# Crear un lollipop plot con ggplot2
g1<- ggplot(top_food, aes(y = Country, x = Perc_ForestDI05)) +
  geom_segment(aes(y = Country, yend = Country, x = 0, xend = Perc_ForestDI05), color = "black") +
  geom_point(size = 5, color = "red", fill = alpha("orange", 0.3), alpha = 0.7, shape = 21, stroke = 2) +
  labs(y = "Country",x = "% forest extent with DI≥0.3") + theme_bw()+ xlim(0,100)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
g1
ggsave(filename="Figures/Maps/Lolipoop_Food.pdf",plot=g1, width=9, height=21, units="cm" )

#### PET MAP LOLIPOOP #########

# Cargar el archivo shapefile de países
countries <- st_read("Data/Spatial/Countries/TropicalCountries.shp")
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
countries <- st_transform(countries, mollweide)

# Cargar el raster
pet_map <- rast("Results/Output/Map_Hunting/pet_map.tif")

# Calcular las dimensiones de los píxeles en km
x <- res(pet_map)[1] / 1000
y <- res(pet_map)[2] / 1000

# Calcular el área de cada píxel en km^2
area_km2 <- x * y

# Iterar sobre cada país
result_list <- list()


for (i in 1:nrow(countries)) {
  # Obtener el polígono del país
  country_polygon <- st_geometry(countries[i, ])
  
  # Convertir el polígono a objeto terra
  terra_poly <- vect(st_as_sf(countries[i, ]))
  
  # Aplicar la máscara al raster usando el polígono
  masked_raster <- mask(pet_map, terra_poly)
  
  # Convertir el raster a data.table
  masked_dt <- as.data.table(masked_raster, long = TRUE)
  
  # Calcular el número total de píxeles forestales para el país
  forest_extent <- sum(!is.na(masked_dt$pet_map), na.rm = TRUE)
  
  # Calcular el número de píxeles con DI > 0.3
  pixels_DI_gt_05 <- sum(masked_dt$pet_map >= 0.3, na.rm = TRUE)
  
  # Calcular la columna Forest_Extent y Forest_km2
  country_df <- data.frame(
    Country = countries$COUNTRY[i],
    Forest_Extent = forest_extent,
    Forest_km2 = forest_extent * area_km2
  )
  
  # Calcular la variable DI
  DI <- pixels_DI_gt_05
  
  # Calcular DI_km2
  DI_km2 <- DI * area_km2
  
  # Calcular Perc_ForestDI05
  Perc_ForestDI05 <- (DI_km2/country_df$Forest_km2) * 100
  
  # Agregar las nuevas columnas al data.frame
  country_df <- cbind(country_df, DI= DI, DI_km2 = DI_km2, Perc_ForestDI05 = Perc_ForestDI05)
  
  # Agregar la lista al resultado
  result_list[[i]] <- country_df
}

# Combinar las listas en un solo data frame
result_df <- do.call(rbind, result_list)

# Ordenar el dataframe por DI_km2 en orden descendente
result_df <- result_df %>% arrange(desc(Perc_ForestDI05))
fwrite(result_df, "Results/Output/Map_Hunting/Pet_Defaunation_Area.csv")


result_df02 <- fread("Results/Output/Map_Hunting/Pet_Defaunation_Area.csv")
result_df02 <- result_df02[-which(result_df02$Country=="United States"),]
# Tomar los 35 países con mayores valores de DI_km2
top_pet <- head(result_df02, 35)
top_pet$Country <- factor(top_pet$Country, levels = top_pet$Country[order(top_pet$Perc_ForestDI05)])

# Crear un lollipop plot con ggplot2
g2<- ggplot(top_pet, aes(y = Country, x = Perc_ForestDI05)) +
  geom_segment(aes(y = Country, yend = Country, x = 0, xend = Perc_ForestDI05), color = "black") +
  geom_point(size = 5, color = "red", fill = alpha("orange", 0.3), alpha = 0.7, shape = 21, stroke = 2) +
  labs(y = "Country",x = "% forest extent with DI≥0.3") + theme_bw()+ xlim(0,100)+
  theme(axis.text.x = element_text(angle = 0, hjust = 1))
g2
ggsave(filename="Figures/Maps/Lolipoop_Pet.pdf",plot=g2, width=9, height=21, units="cm" )




#### ESTIMATE EXTENT OF HOTSPOTS 

China <- mask(pet_map,countries[countries$COUNTRY=="China",])
China <- as.data.table(China)
China_DI05 <- sum(China$pet_map >= 0.3, na.rm = TRUE)
China_DI05 <- China_DI05*area_km2
rm(list=c("China", "China_DI05"))
gc()

Indonesia <- mask(pet_map,countries[countries$COUNTRY=="Indonesia",])
Indonesia <- as.data.table(Indonesia)
Indonesia_DI05 <- sum(Indonesia$pet_map >= 0.3, na.rm = TRUE)
Indonesia_DI05 <- Indonesia_DI05*area_km2
rm(list=c("Indonesia", "Indonesia_DI05"))
gc()

Brazil <- mask(pet_map,countries[countries$COUNTRY=="Brazil",])
Brazil <- as.data.table(Brazil)
Brazil_DI05 <- sum(Brazil$pet_map >= 0.3, na.rm = TRUE)
Brazil_DI05 <- Brazil_DI05*area_km2
rm(list=c("Brazil", "Brazil_DI05"))
gc()


PR <- mask(pet_map,countries[countries$COUNTRY=="Puerto Rico (US)",])
PR <- as.data.table(PR)
PR_DI05 <- sum(PR$pet_map >= 0.3, na.rm = TRUE)
PR_DI05 <- PR_DI05*area_km2
rm(list=c("PR", "PR_DI05"))

Jamaica <- mask(pet_map,countries[countries$COUNTRY=="Jamaica",])
Jamaica <- as.data.table(Jamaica)
Jamaica_DI05 <- sum(Jamaica$pet_map >= 0.3, na.rm = TRUE)
Jamaica_DI05 <- Jamaica_DI05*area_km2
rm(list=c("Jamaica", "Jamaica_DI05"))

Taiwan <- mask(pet_map,countries[countries$COUNTRY=="Taiwan",])
Taiwan <- as.data.table(Taiwan)
Taiwan_DI05 <- sum(Taiwan$pet_map >= 0.3, na.rm = TRUE)
Taiwan_DI05 <- Taiwan_DI05*area_km2
rm(list=c("Taiwan", "Taiwan_DI05"))

TT <- mask(pet_map,countries[countries$COUNTRY=="Trinidad and Tobago",])
TT <- as.data.table(TT)
TT_DI05 <- sum(TT$pet_map >= 0.3, na.rm = TRUE)
TT_DI05 <- TT_DI05*area_km2
rm(list=c("TT", "TT_DI05"))

Gambia <- mask(pet_map,countries[countries$COUNTRY=="Gambia",])
Gambia <- as.data.table(Gambia)
Gambia_DI05 <- sum(Gambia$pet_map >= 0.3, na.rm = TRUE)
Gambia_DI05 <- Gambia_DI05*area_km2
rm(list=c("Gambia", "Gambia_DI05"))

Dom_Rep  <- mask(pet_map,countries[countries$COUNTRY=="Dominican Republic",])
Dom_Rep <- as.data.table(Dom_Rep)
Dom_Rep_DI05 <- sum(Dom_Rep$pet_map >= 0.3, na.rm = TRUE)
Dom_Rep_DI05 <- Dom_Rep_DI05*area_km2
rm(list=c("Dom_Rep", "Dom_Rep_DI05"))


pet_map <- as.data.table(pet_map)

# Sumatra island stat
sumatra <- st_read("Data/Spatial/Countries/Sumatra_Island.shp")
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
sumatra <- st_transform(sumatra, mollweide)
sumatra <- mask(pet_map,sumatra)
sumatra <- as.data.table(sumatra)
median(sumatra$pet_map, na.rm = TRUE)
sd(sumatra$pet_map, na.rm = TRUE)
China_DI05 <- China_DI05*area_km2
rm(list=c("China", "China_DI05"))

# Java island stat
java <- st_read("Data/Spatial/Countries/Java_Island.shp")
mollweide <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
java  <- st_transform(java , mollweide)
java  <- mask(pet_map,java )
java  <- as.data.table(java)
mean(java $pet_map, na.rm = TRUE)
sd(java $pet_map, na.rm = TRUE)
median(java $pet_map, na.rm = TRUE)
java_DI08 <- sum(java$pet_map >= 0.8, na.rm = TRUE)
java_DI08 <- java_DI08*area_km2
java_DI08/nrow(java)*area_km2 
gc()

###### Supplementary table ######
table <- left_join(result_df[1:85,c("Country","Forest_km2", "DI_km2", "Perc_ForestDI05")],
                   result_df02[1:85,c("Country","Forest_km2", "DI_km2", "Perc_ForestDI05")], 
                   by="Country")
table <- arrange(table, Country)
table$Perc_ForestDI05.x <- round(table$Perc_ForestDI05.x,2)
table$Perc_ForestDI05.y <- round(table$Perc_ForestDI05.y,2)
table$Forest_km2.x <- round(table$Forest_km2.x,2)
table$DI_km2.x <- round(table$DI_km2.x,2)
table$DI_km2.y <- round(table$DI_km2.y,2)
writexl::write_xlsx(table,"Results/Output/Map_Hunting/Countries_Pet_Food.xlsx")

 