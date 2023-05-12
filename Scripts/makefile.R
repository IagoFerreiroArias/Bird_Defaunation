#This is a master script in order to run analysis

#This project contains several scripts that should be run in some order

#1: PCA_Birds.R
# Explore diet categories for bird and PCA's axis as potential predictors
source("Script/PCA_Birds.R")

#2: Phylo_Birds.R
# Explore evolutionary distinctiveness as a proxy of bird attractiveness for hunters 
# see Sheffers et al 2019, prob of being traded ~ evo distinctiviness
source("Script/Phylo_Birds.R")

#3: Traits_Birds.R
# Data cleaning of AVONET and Sheffers data 
source("Script/Trait_Birds.R")

#4: Data_Cleaning.R
# Clean database on hunting impacts and merge with previous data frames created (e.g. traits)
source("Script/Data_Cleaning.R")
  
#5:: RR_BRMS
# Modelling hunting impacts using bayessian regression models
source("Script/RR_BRMS.R")

#6: Model_Visualization.R
# Diagnosis plots and marginal effects from hurdle model
source("Script/Model_Visualization.R")

#6: Variable_Importance.R
# Estimation of variable importance from hurdle model.
source("Script/Varibale_Importance.R")





