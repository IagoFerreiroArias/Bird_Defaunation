# Script: Elton traits and PCA diet score of birds. Fixed Effect 
# Author: Iago Ferreiro Arias
# Date: 20th September 2022

###############################################################################
library(dplyr)
library(FactoMineR) #to carry out PCA
library(factoextra) #to extract data from PCA 
library(writexl) #export results into xlsx files
##############################################################################

#Import Elton Traits database
Elton <- read.delim("Data/Traits/ELTON/Elton_BirdFuncDat.txt")
table(Elton$Diet.Certainty)

Elton$Species <- Elton$Scientific
Elton <- Elton[1:9993,] #remove empty rows (2 last rows)

# create matrix for PCA

rownames(Elton) <- make.names(Elton$Species,unique=T) #rownames have points between genus and species specific name
diet_matrix <- Elton %>% dplyr::select(Diet.Inv, Diet.Vend, Diet.Vfish, Diet.Vunk, 
                                       Diet.Scav, Diet.Fruit, Diet.Nect, Diet.Seed,
                                       Diet.PlantO)
diet_matrix<- as.matrix(diet_matrix, rownames.force = T) 

# Run Principal Component Analysis
PCA_diet <- FactoMineR::PCA(diet_matrix, graph=F) # 2 missing values

# Extract coordinates values for each species from PCA
PCA_score <-  as_tibble(factoextra::get_pca_ind(PCA_diet)$coord) 

# Plot PCA results
fviz_pca_biplot(PCA_diet, repel=T, label="var", ggtheme = theme_minimal(),col.ind = "contrib", 
                gradient.cols = c("#FBEC82", "#BFD37D", "#CDFB82","#64B73A", "#A57DD3","#865EC9"))

# Save figure
ggsave(fviz_pca_biplot(PCA_diet, repel=T, label="var", ggtheme = theme_minimal(),col.ind = "contrib", 
                       gradient.cols = c("#FBEC82", "#BFD37D", "#CDFB82","#64B73A", "#A57DD3","#865EC9" )), 
                       path="Figures/PCA Diet", filename="PCA Bird Diet Biplot.pdf")

#PCA coords may be influenced by imputed values of diet
#Check distribution of certainty in diet values across spcies
Elton$Diet.Certainty <- as.factor(Elton$Diet.Certainty)
table(Elton$Diet.Certainty)

# A – Highly certain that the source is reliable and has provided accurate diet information.
# B – Reasonably confident that source is reliable, however information on relative proportion 
# of each diet category making up the whole diet is possibly uncertain.
# C – Quality of source diet estimate unclear or inferred from a specific congeneric species with good sources.
# D1 – Species-level information missing and value is that typical for genus (D1).
# D2 – Species-level information missing and value is that typical for family (D2).


# Contribution of each diet item to PCA dimensions
fviz_contrib(PCA_diet, choice="var", axes= 1) #Axis 1
fviz_contrib(PCA_diet, choice="var", axes= 2) #Axis 2

# Save figure
ggsave(fviz_contrib(PCA_diet, choice="var", axes= 1), path="Figures/PCA Diet", filename="Diet Item Contributions to Dim1.pdf")
ggsave(fviz_contrib(PCA_diet, choice="var", axes= 2), path="Figures/PCA Diet", filename="Diet Item Contributions to Dim2.pdf") 


# Top 10 species contributions to PCA dimensions
fviz_contrib(PCA_diet, choice="ind", axes = 1, top=10) #Axis 1
fviz_contrib(PCA_diet, choice="ind", axes = 2, top=10) #Axis 2

# Save figures
ggsave(fviz_contrib(PCA_diet, choice="ind", axes = 1, top=10), path="Figures/PCA Diet", filename="Top10 Species' Contributions to Dim1.pdf")
ggsave(fviz_contrib(PCA_diet, choice="ind", axes = 2, top=10), path="Figures/PCA Diet", filename="Top10 Species' Contributions to Dim2.pdf")


#Merge Elton database and PCA coordinates
Elton <- read.delim("Data/Traits/ELTON/Elton_BirdFuncDat.txt")
PCA_score <-  as_tibble(factoextra::get_pca_ind(PCA_diet)$coord)
Elton <- Elton[1:9993,] #remove empty rows (2 last rows)
Elton <- cbind(Elton, PCA_score)

#Clean environment
rm(list=c("diet_matrix", "PCA_diet", "PCA_score"))

# Reduce ELTON TRAITS dataset to only predictors that we are interested

Elton$Species <- Elton$Scientific
Elton$PCA_Dim1 <- Elton$Dim.1
Elton$PCA_Dim2 <- Elton$Dim.2
Elton$Activity <- ifelse(Elton$Nocturnal==0, "Diurnal", "Nocturnal")

Elton <- Elton %>% dplyr::select(Species, Activity, PCA_Dim1, PCA_Dim2)

