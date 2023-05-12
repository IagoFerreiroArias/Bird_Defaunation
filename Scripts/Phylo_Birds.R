# Script: Phylogenetic signal on hunting impacts on tropical mammals
# Author: Iago Ferreiro Arias
# Date:  20th September 2022

################################################################################
library(dplyr)
library(writexl) #export results into xlsx files
library(tictoc) 
library(ape) #to extract covariance matrix from phylo trees
library(caper) # phylogenetic signal in binary trait
library(picante) #evolutionary distinctiveness
################################################################################


# Import and prepare consensus tree from Jetz et al 2012. 
# Available for download at http://vertlife.org/phylosubsets

consensus_tree <- ape::read.nexus("Data/Phylo/ConsensusTree150_05credibility.tree")

#Calculate evolutionary distinctiveness for each species
evo_distinctiveness <- picante::evol.distinct(consensus_tree)
evo_distinctiveness <- evo_distinctiveness %>% rename(Evo_Distinctiveness = w)
evo_distinctiveness$Species <- gsub("_", " ", evo_distinctiveness$Species)

#Calculate taxonomic distinctiveness for each species
taxon_dist <- picante::tax.distinctiveness(consensus_tree, type=c("Vane-Wright"))
taxon_dist <- taxon_dist %>% rename(Tax_Distinctiveness = w)
taxon_dist$Species <- gsub("_", " ", taxon_dist$Species)

phylo_traits <- left_join(evo_distinctiveness, taxon_dist, by="Species")

rm(list=c("evo_distinctiveness", "taxon_dist"))

# create phylogenetic covariance matrix for all bird species

bird_cov_matrix <- ape::vcv.phylo(consensus_tree) 
phylo_cov <- bird_cov_matrix 
write.csv(bird_cov_matrix, "Data/Phylo/Jetz_Bird_PhyloCov_Matrix.csv")
rm(bird_cov_matrix)

#####

