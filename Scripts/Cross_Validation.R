# Script: Cross-validation of hurdle model
# Author: Iago Ferreiro- Arias
# Date: 12th May 2023

library(brms) #bayesian hurdle modelling
library(loo) # k-fold crossvalidation
library(ape) #vcov matrix for hurdle models
library(tictoc)
library(phytools) #phylogenetic autocorrelation
library(raster) #create spatial blocks
library(DHARMa) #Spatial autocorrelation
library(caret)
library(dplyr)


RR_data <- read.csv("Data/Bird_RR_data.csv")

# Categorical predictor: non-hunted species vs species hunted (for food and/or trade)
RR_data$Hunted <- ifelse(RR_data$Food == "No" & RR_data$Traded == "No", "No", "Yes")
#unique(RR_data[RR_data$Hunted=="No","BirdTree_Species"])

RR_data[which(RR_data$BirdTree_Species=="Crypturellus bartletti"), "Hunted"] <- "Yes" #scabin
RR_data[which(RR_data$BirdTree_Species=="Ortalis guttata"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Tockus camurus"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Tockus hartlaubi"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Bycanistes albotibialis"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Geotrygon veraguensis"), "Hunted"] <- "Yes" #red list of bird ecuador
RR_data[which(RR_data$BirdTree_Species=="Geotrygon chiriquensis"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Geotrygon costaricensis"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Geotrygon lawrencii"), "Hunted"] <- "Yes"
RR_data[which(RR_data$BirdTree_Species=="Patagioenas nigrirostris"), "Hunted"] <- "Yes" 
RR_data[which(RR_data$BirdTree_Species=="Patagioenas plumbea"), "Hunted"] <- "Yes"

RR_data <- subset(RR_data, Hunted=="Yes")

# Body mass
RR_data$BodyMass_log <- log10(RR_data$Body_Mass/1000)
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)
RR_data$BodyMass_log <- as.numeric(RR_data$BodyMass_log)

# Distance to hunter's access point

RR_data$DistHunt_log <- log10(RR_data$Dist_Hunters +0.1) 
RR_data$DistHunt_log <- scale(RR_data$DistHunt_log, center = T, scale = T)
RR_data$DistHunt_log <- as.numeric(RR_data$DistHunt_log)

# Travel distance to nearest town
RR_data$TravDist_log <- log10(RR_data$TravDist +0.1) 
RR_data$TravDist_log <- scale(RR_data$TravDist_log, center = T, scale = T)
RR_data$TravDist_log <- as.numeric(RR_data$TravDist_log)

# Stunting prevalence
RR_data$Stunting_log <- log10(RR_data$Stunting + 0.1)
RR_data$Stunting_log <- scale(RR_data$Stunting_log, center=T, scale=T)
RR_data$Stunting_log <- as.numeric(RR_data$Stunting_log)

# Human population density
RR_data$PopDens_log <- log10(RR_data$PopDens + 0.1)
RR_data$PopDens_log <- scale(RR_data$PopDens_log, center=T, scale=T)
RR_data$PopDens_log <- as.numeric(RR_data$PopDens_log)

# NPP
RR_data$NPP_log <- log10(RR_data$NPP + 0.1)
RR_data$NPP_log <- scale(RR_data$NPP_log, center=T, scale=T)
RR_data$NPP_log <- as.numeric(RR_data$NPP_log)

# Establish reference levels for categorical variables:
RR_data$Reserve<-as.factor(RR_data$Reserve)
RR_data$Reserve<-relevel(RR_data$Reserve, ref ="0")

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)

RR_data <- subset(RR_data, RR<5) # remove outliers

# Import and prepare consensus tree from Jetz et al 2012. 

consensus_tree <- ape::read.nexus("Data/Phylo/ConsensusTree150_05credibility.tree")

#Create phylogenetic correlation matrix
RR_data$BirdTree_Species <- sub(" ", "_", RR_data$BirdTree_Species)
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label), ]

#Cyanoderma chrysaeum = Stachyris chrysaea
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label),"BirdTree_Species"] <- "Stachyris_chrysaea"

drops <-consensus_tree$tip.label[!consensus_tree$tip.label %in% RR_data$BirdTree_Species]
phylo_tree <- drop.tip(consensus_tree, drops)
phylo_cov <- ape::vcv.phylo(phylo_tree, corr=FALSE) #correlation matrix. If corr= TRUE => vcv matrix
#rm(list=c("drops", "phylo_tree", "consensus_tree"))

pr <-c(prior(normal(0, 10), class='Intercept'), prior(normal(0, 1), class='b'))#weakly informative priors


hurdle_formula <- bf(RR ~ DistHunt_log + BodyMass_log + Stunting_log + TravDist_log + 
                       PopDens_log +  NPP_log + Reserve + 
                       BodyMass_log:DistHunt_log + BodyMass_log:TravDist_log + 
                       BodyMass_log:Stunting_log + NPP_log:DistHunt_log +
                       (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|CountryNum),
                     hu ~ DistHunt_log + BodyMass_log + Stunting_log + TravDist_log + 
                       PopDens_log +  NPP_log + Reserve + 
                       BodyMass_log:DistHunt_log + BodyMass_log:TravDist_log + 
                       BodyMass_log:Stunting_log + NPP_log:DistHunt_log  +
                       (1|CountryNum) + (1|gr(BirdTree_Species, cov = phylo_cov)))

set.seed(123) #reproducible analysis

RR_data$Def <- ifelse(RR_data$RR<0.3, "High",
                      ifelse(RR_data$RR>0.3 & RR_data$RR <0.7, "Moderate",
                             ifelse(RR_data$RR>0.7 & RR_data$RR <1, "Low", "Increase")))

####### Phylogenetic cross-validation: Orders  ########

# When there are no positive results, sensitivity and specificity is not defined and a value of NA is returned.
# by blocking by orders, families, species or spatially, some blocks may not have some RR in those categories and
# some NAs are returned

RMSE<-rep(NA,10)
SensHigh<-rep(NA,10)
SensMod <-rep(NA,10)
SensLow <-rep(NA,10)
SensInc<-rep(NA,10)

SpecHigh<-rep(NA,10)
SpecMod <-rep(NA,10)
SpecLow<-rep(NA,10)
SpecInc<-rep(NA,10)

AccHigh<-rep(NA,10)
AccMod<-rep(NA,10)
AccLow<-rep(NA,10)
AccInc<-rep(NA,10)

RR_data$foldsOrder <- kfold_split_grouped(K=10, x=RR_data$Order)

tic()
for (i in unique(RR_data$foldsOrder)){
  
  testing<-RR_data[RR_data$foldsOrder %in% i,]
  training<-RR_data[!RR_data$foldsOrder %in% i,]
  Reserve<-all(unique(testing$Reserve) %in% unique(training$Reserve)) #always true?
  if(!all(Reserve)) {next}
  
  mod <- brm(hurdle_formula,
             data = RR_data,  
             chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
             family=hurdle_lognormal(), prior=pr,  
             data2 = list(phylo_cov = phylo_cov), seed=123,
             control = list(adapt_delta =0.99, max_treedepth=15))
  
  pred<-predict(mod, newdata=testing, summary=TRUE, allow_new_levels=TRUE, type= "response")
  RMSE[i]<-sqrt(mean((pred[,1] - testing$RR)^2)) # Root Mean Square Erorr
  
  confusion_df <- data.frame("Def"=testing$Def, "Pred"=pred[,1])
  confusion_df$Def_pred <- ifelse(confusion_df$Pred<0.3, "High",
                           ifelse(confusion_df$Pred>0.3 & confusion_df$Pred <0.7, "Moderate",
                           ifelse(confusion_df$Pred>0.7 & confusion_df$Pred <1, "Low", "Increase")))
  
  confusion_df$Def <- factor(confusion_df$Def, levels=c("High", "Moderate", "Low","Increase"))
  confusion_df$Def_pred <- factor(confusion_df$Def_pred, levels=c("High", "Moderate", "Low","Increase"))
  
  # Sensitivity: True Positive / True positive + False Negative
  # Specificity: True Negative / True Negative + False Positive
  
  #caret functions works with two levels factors
  
  #Sensitivity and Specificity for high impacts
  
  confusion_df$Def_High <- ifelse(confusion_df$Def=="High", "Yes", "No")
  confusion_df$Pred_High <- ifelse(confusion_df$Def_pred=="High", "Yes", "No")
  confusion_df$Def_High <- factor(confusion_df$Def_High, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_High <- factor(confusion_df$Pred_High, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensHigh[i]<- caret::sensitivity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  SpecHigh[i]<- caret::specificity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  #AccHigh <- (SensHigh[i] + SpecHigh[i])/2
  
  #Sensitivity and Specificity for moderate impacts
  
  confusion_df$Def_Mod <- ifelse(confusion_df$Def=="Moderate", "Yes", "No")
  confusion_df$Pred_Mod <- ifelse(confusion_df$Def_pred=="Moderate", "Yes", "No")
  
  confusion_df$Def_Mod <- factor(confusion_df$Def_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Mod <- factor(confusion_df$Pred_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensMod[i]<- caret::sensitivity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  SpecMod[i]<- caret::specificity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  #AccMod[i] <- (SensMod[i] + SpecMod[i])/2
  
  #Sensitivity and Specificity for low impacts
  
  confusion_df$Def_Low <- ifelse(confusion_df$Def=="Low", "Yes", "No")
  confusion_df$Pred_Low <- ifelse(confusion_df$Def_pred=="Low", "Yes", "No")
  
  confusion_df$Def_Low <- factor(confusion_df$Def_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Low <- factor(confusion_df$Pred_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensLow[i]<- caret::sensitivity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  SpecLow[i]<- caret::specificity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  #AccLow[i] <- (SensLow[i] + SpecLow[i])/2
  
  #Sensitivity and Specificity for increase in abundance
  
  confusion_df$Def_Inc <- ifelse(confusion_df$Def=="Increase", "Yes", "No")
  confusion_df$Pred_Inc  <- ifelse(confusion_df$Def_pred=="Increase", "Yes", "No")
  
  confusion_df$Def_Inc <- factor(confusion_df$Def_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Inc <- factor(confusion_df$Pred_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensInc[i]<- caret::sensitivity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  SpecInc[i]<- caret::specificity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  #AccInc[i] <- (SensInc[i] + SpecInc[i])/2
  
}

toc()

Cross_Order <- data.frame("Fold"=seq(1,10,1), "Block" = rep("Order",10))
Cross_Order <- cbind(Cross_Order, RMSE, 
                     SensHigh, SensMod, SensLow, SensInc,
                     SpecHigh, SpecMod, SpecLow, SpecInc, 
                     AccHigh, AccMod, AccLow, AccInc)

#write.csv(Cross_Order, 'Results/Cross_Order.csv', row.names=FALSE)

#### Phylogenetic cross-validation: Family ########

RMSE<-rep(NA,10)
SensHigh<-rep(NA,10)
SensMod <-rep(NA,10)
SensLow <-rep(NA,10)
SensInc<-rep(NA,10)

SpecHigh<-rep(NA,10)
SpecMod <-rep(NA,10)
SpecLow<-rep(NA,10)
SpecInc<-rep(NA,10)

AccHigh<-rep(NA,10)
AccMod<-rep(NA,10)
AccLow<-rep(NA,10)
AccInc<-rep(NA,10)


RR_data$foldsFamily <- kfold_split_grouped(K=10, x=RR_data$Family)

tic()
for (i in unique(RR_data$foldsFamily )){
  
  testing<-RR_data[RR_data$foldsFamily  %in% i,]
  training<-RR_data[!RR_data$foldsFamily  %in% i,]
  Reserve<-all(unique(testing$Reserve) %in% unique(training$Reserve)) #always true?
  if(!all(Reserve)) {next}
  
  mod <- brm(hurdle_formula,
             data = RR_data,  
             chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
             family=hurdle_lognormal(), prior=pr,  
             data2 = list(phylo_cov = phylo_cov), seed=123,
             control = list(adapt_delta =0.99, max_treedepth=15))
  
  pred<-predict(mod, newdata=testing, summary=TRUE, allow_new_levels=TRUE, type= "response")
  RMSE[i]<-sqrt(mean((pred[,1] - testing$RR)^2)) # Root Mean Square Error
  
  confusion_df <- data.frame("Def"=testing$Def, "Pred"=pred[,1])
  confusion_df$Def_pred <- ifelse(confusion_df$Pred<0.3, "High",
                                  ifelse(confusion_df$Pred>0.3 & confusion_df$Pred <0.7, "Moderate",
                                         ifelse(confusion_df$Pred>0.7 & confusion_df$Pred <1, "Low", "Increase")))
  
  confusion_df$Def <- factor(confusion_df$Def, levels=c("High", "Moderate", "Low","Increase"))
  confusion_df$Def_pred <- factor(confusion_df$Def_pred, levels=c("High", "Moderate", "Low","Increase"))
  
  # Sensitivity: True Positive / True positive + False Negative
  # Specificity: True Negative / True Negative + False Positive
  
  #caret functions works with two levels factors
  
  #Sensitivity and Specificity for high impacts
  
  confusion_df$Def_High <- ifelse(confusion_df$Def=="High", "Yes", "No")
  confusion_df$Pred_High <- ifelse(confusion_df$Def_pred=="High", "Yes", "No")
  confusion_df$Def_High <- factor(confusion_df$Def_High, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_High <- factor(confusion_df$Pred_High, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensHigh[i]<- caret::sensitivity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  SpecHigh[i]<- caret::specificity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  AccHigh[i] <- (SensHigh[i] + SpecHigh[i])/2
  
  #Sensitivity and Specificity for moderate impacts
  
  confusion_df$Def_Mod <- ifelse(confusion_df$Def=="Moderate", "Yes", "No")
  confusion_df$Pred_Mod <- ifelse(confusion_df$Def_pred=="Moderate", "Yes", "No")
  
  confusion_df$Def_Mod <- factor(confusion_df$Def_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Mod <- factor(confusion_df$Pred_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensMod[i]<- caret::sensitivity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  SpecMod[i]<- caret::specificity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  AccMod[i] <- (SensMod[i] + SpecMod[i])/2
  
  #Sensitivity and Specificity for low impacts
  
  confusion_df$Def_Low <- ifelse(confusion_df$Def=="Low", "Yes", "No")
  confusion_df$Pred_Low <- ifelse(confusion_df$Def_pred=="Low", "Yes", "No")
  
  confusion_df$Def_Low <- factor(confusion_df$Def_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Low <- factor(confusion_df$Pred_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensLow[i]<- caret::sensitivity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  SpecLow[i]<- caret::specificity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  AccLow[i] <- (SensLow[i] + SpecLow[i])/2
  
  #Sensitivity and Specificity for increase in abundance
  
  confusion_df$Def_Inc <- ifelse(confusion_df$Def=="Increase", "Yes", "No")
  confusion_df$Pred_Inc  <- ifelse(confusion_df$Def_pred=="Increase", "Yes", "No")
  
  confusion_df$Def_Inc <- factor(confusion_df$Def_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Inc <- factor(confusion_df$Pred_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensInc[i]<- caret::sensitivity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  SpecInc[i]<- caret::specificity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  AccInc[i] <- (SensInc[i] + SpecInc[i])/2

}

Cross_Family <- data.frame("Fold"=seq(1,10,1), "Block" = rep("Family",10))
Cross_Family <- cbind(Cross_Family, RMSE, 
                      SensHigh, SensMod, SensLow, SensInc,
                      SpecHigh, SpecMod, SpecLow, SpecInc, 
                      AccHigh, AccMod, AccLow, AccInc)

#write.csv(Cross_Family, 'Results/Cross_Family.csv', row.names=FALSE)

######## Phylogenetic cross-validation: Species ##########

RMSE<-rep(NA,10)
SensHigh<-rep(NA,10)
SensMod <-rep(NA,10)
SensLow <-rep(NA,10)
SensInc<-rep(NA,10)

SpecHigh<-rep(NA,10)
SpecMod <-rep(NA,10)
SpecLow<-rep(NA,10)
SpecInc<-rep(NA,10)

AccHigh<-rep(NA,10)
AccMod<-rep(NA,10)
AccLow<-rep(NA,10)
AccInc<-rep(NA,10)

RR_data$foldsSpecies <- kfold_split_grouped(K=10, x=RR_data$BirdTree_Species)

tic()
for (i in unique(RR_data$foldsSpecies)){
  
  testing<-RR_data[RR_data$foldsSpecies  %in% i,]
  training<-RR_data[!RR_data$foldsSpecies %in% i,]
  Reserve<-all(unique(testing$Reserve) %in% unique(training$Reserve)) #always true?
  if(!all(Reserve)) {next}
  
  mod <- brm(hurdle_formula,
             data = RR_data,  
             chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
             family=hurdle_lognormal(), prior=pr,  
             data2 = list(phylo_cov = phylo_cov), seed=123,
             control = list(adapt_delta =0.99, max_treedepth=15))
  
  pred<-predict(mod, newdata=testing, summary=TRUE, allow_new_levels=TRUE, type= "response")
  RMSE[i]<-sqrt(mean((pred[,1] - testing$RR)^2)) # Root Mean Square Error
  
  confusion_df <- data.frame("Def"=testing$Def, "Pred"=pred[,1])
  confusion_df$Def_pred <- ifelse(confusion_df$Pred<0.3, "High",
                                  ifelse(confusion_df$Pred>0.3 & confusion_df$Pred <0.7, "Moderate",
                                         ifelse(confusion_df$Pred>0.7 & confusion_df$Pred <1, "Low", "Increase")))
  
  confusion_df$Def <- factor(confusion_df$Def, levels=c("High", "Moderate", "Low","Increase"))
  confusion_df$Def_pred <- factor(confusion_df$Def_pred, levels=c("High", "Moderate", "Low","Increase"))
  
  # Sensitivity: True Positive / True positive + False Negative
  # Specificity: True Negative / True Negative + False Positive
  
  #caret functions works with two levels factors
  
  #Sensitivity and Specificity for high impacts
  
  confusion_df$Def_High <- ifelse(confusion_df$Def=="High", "Yes", "No")
  confusion_df$Pred_High <- ifelse(confusion_df$Def_pred=="High", "Yes", "No")
  confusion_df$Def_High <- factor(confusion_df$Def_High, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_High <- factor(confusion_df$Pred_High, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensHigh[i]<- caret::sensitivity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  SpecHigh[i]<- caret::specificity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  AccHigh[i] <- (SensHigh[i] + SpecHigh[i])/2
  
  #Sensitivity and Specificity for moderate impacts
  
  confusion_df$Def_Mod <- ifelse(confusion_df$Def=="Moderate", "Yes", "No")
  confusion_df$Pred_Mod <- ifelse(confusion_df$Def_pred=="Moderate", "Yes", "No")
  
  confusion_df$Def_Mod <- factor(confusion_df$Def_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Mod <- factor(confusion_df$Pred_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensMod[i]<- caret::sensitivity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  SpecMod[i]<- caret::specificity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  AccMod[i] <- (SensMod[i] + SpecMod[i])/2
  
  #Sensitivity and Specificity for low impacts
  
  confusion_df$Def_Low <- ifelse(confusion_df$Def=="Low", "Yes", "No")
  confusion_df$Pred_Low <- ifelse(confusion_df$Def_pred=="Low", "Yes", "No")
  
  confusion_df$Def_Low <- factor(confusion_df$Def_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Low <- factor(confusion_df$Pred_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensLow[i]<- caret::sensitivity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  SpecLow[i]<- caret::specificity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  AccLow[i] <- (SensLow[i] + SpecLow[i])/2
  
  #Sensitivity and Specificity for increase in abundance
  
  confusion_df$Def_Inc <- ifelse(confusion_df$Def=="Increase", "Yes", "No")
  confusion_df$Pred_Inc  <- ifelse(confusion_df$Def_pred=="Increase", "Yes", "No")
  
  confusion_df$Def_Inc <- factor(confusion_df$Def_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Inc <- factor(confusion_df$Pred_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensInc[i]<- caret::sensitivity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  SpecInc[i]<- caret::specificity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  AccInc[i] <- (SensInc[i] + SpecInc[i])/2
  
  
}

Cross_Species <- data.frame("Fold"=seq(1,10,1), "Block" = rep("Species",10))
Cross_Species <- cbind(Cross_Species, RMSE, 
                       SensHigh, SensMod, SensLow, SensInc,
                       SpecHigh, SpecMod, SpecLow, SpecInc, 
                       AccHigh, AccMod, AccLow, AccInc)
#write.csv(Cross_Species, 'Results/Cross_Species.csv', row.names=FALSE)

############# Phylogenetic autocorrelation #####################

load("Results/hurdle_model.Rdata")
resid<-residuals(hurdle_model)
residAgg<-tapply(resid[,'Estimate'], RR_data$BirdTree_Species, mean)
Lambda<-phylosig(phylo_tree, residAgg, method='lambda', test=TRUE)
plot(Lambda)
Lambda<-data.frame(Lambda=Lambda$lambda, LogLikelihood=Lambda$logL0, Pvalue=Lambda$P)
write.csv(Lambda,"Results/PhyloSignal_Lambda.csv", row.names=FALSE)

rm(list=c("resid", "residAgg", "Lambda", "phylo_tree", "testing", "training", "phylo_cov")) #clean environment

########### Spatial-block cross-validation for hurdle model ############

blockSize <- 5 #degrees
Xvec<-seq(-180, 180, by=blockSize)
Yvec<-seq(-38, 43, by=blockSize) #Degrees of latitude forest extension: - 37.44 , 42.60 (gis layer)
grid<-expand.grid(Longitude=Xvec, Latitude=Yvec)
grid$ID<-1:nrow(grid)
SpatialBlocks<-raster::rasterFromXYZ(grid)
RR_data$SpatialBlock<-raster::extract(SpatialBlocks, RR_data[,c('Longitude','Latitude')])

RMSE<-rep(NA,10)
SensHigh<-rep(NA,10)
SensMod <-rep(NA,10)
SensLow <-rep(NA,10)
SensInc<-rep(NA,10)

SpecHigh<-rep(NA,10)
SpecMod <-rep(NA,10)
SpecLow<-rep(NA,10)
SpecInc<-rep(NA,10)

AccHigh<-rep(NA,10)
AccMod<-rep(NA,10)
AccLow<-rep(NA,10)
AccInc<-rep(NA,10)

RR_data$foldsSpatial <- kfold_split_grouped(K=10, x=RR_data$SpatialBlock)

for (i in unique(RR_data$foldsSpatial)){
  
  testing<-RR_data[RR_data$foldsSpatial %in% i,]
  training<-RR_data[!RR_data$foldsSpatial %in% i,]
  Reserve<-all(unique(testing$Reserve) %in% unique(training$Reserve)) #always true?
  if(!all(Reserve)) {next}
  
  mod <- brm(hurdle_formula,
             data = RR_data,  
             chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
             family=hurdle_lognormal(), prior=pr,  
             data2 = list(phylo_cov = phylo_cov), seed=123,
             control = list(adapt_delta =0.99, max_treedepth=15))
  
  pred<-predict(mod, newdata=testing, summary=TRUE)
  RMSE[i]<-sqrt(mean((pred[,1] - testing$RR)^2)) # Root Mean Square Error
  
  confusion_df <- data.frame("Def"=testing$Def, "Pred"=pred[,1])
  confusion_df$Def_pred <- ifelse(confusion_df$Pred<0.3, "High",
                                  ifelse(confusion_df$Pred>0.3 & confusion_df$Pred <0.7, "Moderate",
                                         ifelse(confusion_df$Pred>0.7 & confusion_df$Pred <1, "Low", "Increase")))
  
  confusion_df$Def <- factor(confusion_df$Def, levels=c("High", "Moderate", "Low","Increase"))
  confusion_df$Def_pred <- factor(confusion_df$Def_pred, levels=c("High", "Moderate", "Low","Increase"))
  
  # Sensitivity: True Positive / True positive + False Negative
  # Specificity: True Negative / True Negative + False Positive
  
  #caret functions works with two levels factors
  
  #Sensitivity and Specificity for high impacts
  
  confusion_df$Def_High <- ifelse(confusion_df$Def=="High", "Yes", "No")
  confusion_df$Pred_High <- ifelse(confusion_df$Def_pred=="High", "Yes", "No")
  confusion_df$Def_High <- factor(confusion_df$Def_High, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_High <- factor(confusion_df$Pred_High, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensHigh[i]<- caret::sensitivity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  SpecHigh[i]<- caret::specificity(confusion_df$Def_High, confusion_df$Pred_High, positive = levels(confusion_df$Def_High)[1])
  AccHigh[i] <- (SensHigh[i] + SpecHigh[i])/2
  
  #Sensitivity and Specificity for moderate impacts
  
  confusion_df$Def_Mod <- ifelse(confusion_df$Def=="Moderate", "Yes", "No")
  confusion_df$Pred_Mod <- ifelse(confusion_df$Def_pred=="Moderate", "Yes", "No")
  
  confusion_df$Def_Mod <- factor(confusion_df$Def_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Mod <- factor(confusion_df$Pred_Mod, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensMod[i]<- caret::sensitivity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  SpecMod[i]<- caret::specificity(confusion_df$Def_Mod, confusion_df$Pred_Mod, positive = levels(confusion_df$Def_Mod)[1])
  AccMod[i] <- (SensMod[i] + SpecMod[i])/2
  
  #Sensitivity and Specificity for low impacts
  
  confusion_df$Def_Low <- ifelse(confusion_df$Def=="Low", "Yes", "No")
  confusion_df$Pred_Low <- ifelse(confusion_df$Def_pred=="Low", "Yes", "No")
  
  confusion_df$Def_Low <- factor(confusion_df$Def_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Low <- factor(confusion_df$Pred_Low, levels=c("Yes", "No")) #set order of levels for caret's function
  
  SensLow[i]<- caret::sensitivity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  SpecLow[i]<- caret::specificity(confusion_df$Def_Low, confusion_df$Pred_Low, positive = levels(confusion_df$Def_Low)[1])
  AccLow[i] <- (SensLow[i] + SpecLow[i])/2
  
  #Sensitivity and Specificity for increase in abundance
  
  confusion_df$Def_Inc <- ifelse(confusion_df$Def=="Increase", "Yes", "No")
  confusion_df$Pred_Inc  <- ifelse(confusion_df$Def_pred=="Increase", "Yes", "No")
  
  confusion_df$Def_Inc <- factor(confusion_df$Def_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  confusion_df$Pred_Inc <- factor(confusion_df$Pred_Inc, levels=c("Yes", "No")) #set order of levels for caret's function
  #http://127.0.0.1:14263/graphics/plot_zoom_png?width=762&height=832
  SensInc[i]<- caret::sensitivity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  SpecInc[i]<- caret::specificity(confusion_df$Def_Inc, confusion_df$Pred_Inc, positive = levels(confusion_df$Def_Inc)[1])
  AccInc[i] <- (SensInc[i] + SpecInc[i])/2
  
}

Cross_Spatial <- data.frame("Fold"=seq(1,10,1), "Block" = rep("Spatial",10))
Cross_Spatial <- cbind(Cross_Spatial, RMSE, 
                       SensHigh, SensMod, SensLow, SensInc,
                       SpecHigh, SpecMod, SpecLow, SpecInc, 
                       AccHigh, AccMod, AccLow, AccInc)

#write.csv(Cross_Spatial, 'Results/Cross_Spatial.csv', row.names=FALSE)

cross <- rbind(Cross_Order, Cross_Family, Cross_Species,Cross_Spatial)
write.csv(cross, 'Results/CrossBlocks.csv', row.names=FALSE)


cross <- read.csv("Results/CrossBlocks.csv")

RMSE <- cross %>% group_by(Block) %>% summarise(
  RMSE_Mean = mean(RMSE, na.rm=TRUE),
  RMSE_sd = sd(RMSE, na.rm=TRUE)) 

RMSE$Block <- factor(RMSE$Block, levels=c("Order","Family","Species","Spatial"))

ggplot(RMSE) +
     geom_bar( aes(x=Block, y=RMSE_Mean), stat="identity", fill="#4A639B", alpha=0.9) +
     geom_errorbar( aes(x=Block, ymin=RMSE_Mean-RMSE_sd, ymax=RMSE_Mean+RMSE_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
     coord_flip() + theme_classic() + xlab("Block") + ylab("RMSE")

Sensitivity <- cross %>% group_by(Block) %>% summarise(
               SensHigh_Mean = mean(SensHigh, na.rm=TRUE),
               SensMod_Mean = mean(SensMod, na.rm=TRUE),
               SensLow_Mean = mean(SensLow, na.rm=TRUE),
               SensInc_Mean = mean(SensInc, na.rm=TRUE),
               SensHigh_sd = sd(SensHigh, na.rm=TRUE),
               SensMod_sd = sd(SensMod, na.rm=TRUE),
               SensLow_sd = sd(SensLow, na.rm=TRUE),
               SensInc_sd = sd(SensInc, na.rm=TRUE))

Sensitivity$Block <- factor(Sensitivity$Block, levels=c("Order","Family","Species","Spatial"))

s1 <- ggplot(Sensitivity) +
      geom_bar( aes(x=Block, y=SensHigh_Mean), stat="identity", fill="#AEAEAE", alpha=0.9) +
     geom_errorbar( aes(x=Block, ymin=SensHigh_Mean-SensHigh_sd, ymax=SensHigh_Mean+SensHigh_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
     coord_flip() + theme_classic() + xlab("Block") + ylab("Sensitivity") + ylim(0,1.25)+ 
      scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))

s2 <- ggplot(Sensitivity) +
  geom_bar( aes(x=Block, y=SensMod_Mean), stat="identity", fill="#AEAEAE", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SensMod_Mean-SensMod_sd, ymax=SensMod_Mean+SensMod_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Sensitivity") + ylim(0,1)

s3 <- ggplot(Sensitivity) +
  geom_bar( aes(x=Block, y=SensLow_Mean), stat="identity", fill="#AEAEAE", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SensLow_Mean-SensLow_sd, ymax=SensLow_Mean+SensLow_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Sensitivity")+ ylim(-0.1,1)

s4 <- ggplot(Sensitivity) +
  geom_bar( aes(x=Block, y=SensInc_Mean), stat="identity", fill="#AEAEAE", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SensInc_Mean-SensInc_sd, ymax=SensInc_Mean+SensInc_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Sensitivity") + ylim(0,1)


Specificity <- cross %>% group_by(Block) %>% summarise(
  SpecHigh_Mean = mean(SpecHigh, na.rm=TRUE),
  SpecMod_Mean = mean(SpecMod, na.rm=TRUE),
  SpecLow_Mean = mean(SpecLow, na.rm=TRUE),
  SpecInc_Mean = mean(SpecInc, na.rm=TRUE),
  SpecHigh_sd = sd(SpecHigh, na.rm=TRUE),
  SpecMod_sd = sd(SpecMod, na.rm=TRUE),
  SpecLow_sd = sd(SpecLow, na.rm=TRUE),
  SpecInc_sd = sd(SpecInc, na.rm=TRUE))

Specificity$Block <- factor(Specificity$Block, levels=c("Order","Family","Species","Spatial"))

s5 <- ggplot(Specificity) +
  geom_bar( aes(x=Block, y=SpecHigh_Mean), stat="identity", fill="#7B7B7B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SpecHigh_Mean-SpecHigh_sd, ymax=SpecHigh_Mean+SpecHigh_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Specificity") + ylim(0,1)

s6 <- ggplot(Specificity) +
  geom_bar( aes(x=Block, y=SpecMod_Mean), stat="identity", fill="#7B7B7B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SpecMod_Mean-SpecMod_sd, ymax=SpecMod_Mean+SpecMod_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Specificity")+ ylim(0,1)

s7 <- ggplot(Specificity) +
  geom_bar( aes(x=Block, y=SpecLow_Mean), stat="identity", fill="#7B7B7B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SpecLow_Mean-SpecLow_sd, ymax=SpecLow_Mean+SpecLow_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Specificity")+ ylim(0,1)

s8 <- ggplot(Specificity) +
  geom_bar( aes(x=Block, y=SpecInc_Mean), stat="identity", fill="#7B7B7B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=SpecInc_Mean-SpecInc_sd, ymax=SpecInc_Mean+SpecInc_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Specificity")+ ylim(0,1)

ggpubr::ggarrange(s1,s2,s3,s4,s5,s6,s7,s8, ncol=4, nrow=2)


Accuracy <- cross %>% group_by(Block) %>% summarise(
  AccHigh_Mean = mean(AccHigh, na.rm=TRUE),
  AccMod_Mean = mean(AccMod, na.rm=TRUE),
  AccLow_Mean = mean(AccLow, na.rm=TRUE),
  AccInc_Mean = mean(AccInc, na.rm=TRUE),
  AccHigh_sd = sd(AccHigh, na.rm=TRUE),
  AccMod_sd = sd(AccMod, na.rm=TRUE),
  AccLow_sd = sd(AccLow, na.rm=TRUE),
  AccInc_sd = sd(AccInc, na.rm=TRUE))
Accuracy$Block <- factor(Accuracy$Block, levels=c("Order","Family","Species","Spatial"))


Accuracy[Accuracy$Block=="Order","AccHigh_Mean"]<- (Sensitivity[Sensitivity$Block=="Order","SensHigh_Mean"] + 
                                                    Specificity[Specificity$Block=="Order","SpecHigh_Mean"])/2

Accuracy[Accuracy$Block=="Order","AccMod_Mean"]<- (Sensitivity[Sensitivity$Block=="Order","SensMod_Mean"] + 
                                                      Specificity[Specificity$Block=="Order","SpecMod_Mean"])/2

Accuracy[Accuracy$Block=="Order","AccLow_Mean"]<- (Sensitivity[Sensitivity$Block=="Order","SensLow_Mean"] + 
                                                      Specificity[Specificity$Block=="Order","SpecLow_Mean"])/2

Accuracy[Accuracy$Block=="Order","AccInc_Mean"]<- (Sensitivity[Sensitivity$Block=="Order","SensInc_Mean"] + 
                                                      Specificity[Specificity$Block=="Order","SpecInc_Mean"])/2

Accuracy[Accuracy$Block=="Order","AccHigh_sd"]<- (Sensitivity[Sensitivity$Block=="Order","SensHigh_sd"] + 
                                                      Specificity[Specificity$Block=="Order","SpecHigh_sd"])/2

Accuracy[Accuracy$Block=="Order","AccMod_sd"]<- (Sensitivity[Sensitivity$Block=="Order","SensMod_sd"] + 
                                                     Specificity[Specificity$Block=="Order","SpecMod_sd"])/2

Accuracy[Accuracy$Block=="Order","AccLow_sd"]<- (Sensitivity[Sensitivity$Block=="Order","SensLow_sd"] + 
                                                     Specificity[Specificity$Block=="Order","SpecLow_sd"])/2

Accuracy[Accuracy$Block=="Order","AccInc_sd"]<- (Sensitivity[Sensitivity$Block=="Order","SensInc_sd"] + 
                                                     Specificity[Specificity$Block=="Order","SpecInc_sd"])/2

s9 <- ggplot(Accuracy) +
  geom_bar( aes(x=Block, y=AccHigh_Mean), stat="identity", fill="#5B5B5B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=AccHigh_Mean-AccHigh_sd, ymax=AccHigh_Mean+AccHigh_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Balanced Accuracy") + ylim(0,1.1)+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1))

s10 <- ggplot(Accuracy) +
  geom_bar( aes(x=Block, y=AccMod_Mean), stat="identity", fill="#5B5B5B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=AccMod_Mean-AccMod_sd, ymax=AccMod_Mean+AccMod_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Balanced Accuracy")+ ylim(0,1)

s11 <- ggplot(Accuracy) +
  geom_bar( aes(x=Block, y=AccLow_Mean), stat="identity", fill="#5B5B5B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=AccLow_Mean-AccLow_sd, ymax=AccLow_Mean+AccLow_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Balanced Accuracy")+ ylim(0,1)

s12 <- ggplot(Accuracy) +
  geom_bar( aes(x=Block, y=AccInc_Mean), stat="identity", fill="#5B5B5B", alpha=0.9) +
  geom_errorbar( aes(x=Block, ymin=AccInc_Mean-AccInc_sd, ymax=AccInc_Mean+AccInc_sd), linewidth=0.9,width=0.4, colour="orange", alpha=0.9) +
  coord_flip() + theme_classic() + xlab("Block") + ylab("Balanced Accuracy")+ ylim(0,1)


ggpubr::ggarrange(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12, ncol=4, nrow=3)

#clean environment
rm(list=c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10","s11","s12",
          "Sensitivity", "Specificity", "Accuracy", "RMSE"))

# Check spatial autocorrelation model residuals

RR_data$XY<-paste0(RR_data$Longitude,'_',RR_data$Latitude)
XY<-unique(RR_data$XY)
XY2<-matrix(as.numeric(unlist(strsplit(XY, '_'))), ncol=2, byrow=TRUE)
X<-XY2[,1]
Y<-XY2[,2]

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(hurdle_model)),
  observedResponse = RR_data$RR,
  fittedPredictedResponse = apply(t(posterior_epred(hurdle_model)), 1, mean),
  integerResponse = TRUE)

modelcheck_SP<-recalculateResiduals(model_check, group=RR_data$XY)
MoranI<-testSpatialAutocorrelation(modelcheck_SP, x=X, y=Y) #p>0.05 means no spatial autocorrelation
MoranI<-data.frame(Observed=MoranI$statistic[1], Expected=MoranI$statistic[2], SD=MoranI$statistic[3], Pvalue=MoranI$p.value)
row.names(MoranI)<- "MoranI"
write.csv(MoranI,'Results/MoranI.csv', row.names=FALSE)


#rm(list=c("blockSize", "i", "rmse", "X", "Xvec", "Y", "Yvec", "XY", "XY2",
#          "MoranI", "modelcheck_SP", "model_check", "SpatialBlocks", "grid",
#          "Reserve", "Traded"))


