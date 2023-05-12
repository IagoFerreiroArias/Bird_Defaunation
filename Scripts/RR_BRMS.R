<<<<<<< Updated upstream
# Script: Bayesian modelling: Bird's ratio responses to hunting
# Author: Iago Ferreiro Arias
# Date: 3th October 2022

###############################################################################
library(brms)
library(ape)
library(spdep)
library(DHARMa)
library(ggplot2)
library(tictoc)
library(corrplot)
library(dplyr)
library(lme4)
library(performance)
################################################################################

RR_data <- read.csv("Data/Bird_RR_data.csv")

# Log10 and Scale predictors:

# Body mass
RR_data$BodyMass_log <- log10(RR_data$Body_Mass)
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)

#extract  mean and sd to rescale the same variable in the prediction dataset
# Rescaling precitors for projections: data$predictor <- (data$predictor - predictorMean)/predictorSD

df_rescale <- data.frame(Predictor=c("BodyMass", "HWI","DistHunt", "TravDist",
                                     "FoodBiomass", "Stunting", "PopDens"), 
                         Mean=NA, SD=NA)

df_rescale[which(df_rescale$Predictor=="BodyMass") ,"Mean"]<- attr(RR_data$BodyMass_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="BodyMass") ,"SD"]<- attr(RR_data$BodyMass_log, "scaled:scale") #sd

RR_data$BodyMass_log <- as.numeric(RR_data$BodyMass_log)


# Hand Wing Index
RR_data$HWI_log <- log10(RR_data$HWI)
RR_data$HWI_log <- scale(RR_data$HWI_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="HWI") ,"Mean"]<- attr(RR_data$HWI_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="HWI") ,"SD"]<- attr(RR_data$HWI_log, "scaled:scale") #sd

RR_data$HWI_log <- as.numeric(RR_data$HWI_log)


# Distance to hunter's access point

RR_data$DistHunt_log <- log10(RR_data$Dist_Hunters +0.1) 
RR_data$DistHunt_log <- scale(RR_data$DistHunt_log, center = T, scale = T)

df_rescale[which(df_rescale$Predictor=="DistHunt") ,"Mean"]<- attr(RR_data$DistHunt_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="DistHunt") ,"SD"]<- attr(RR_data$DistHunt_log, "scaled:scale") #sd

RR_data$DistHunt_log <- as.numeric(RR_data$DistHunt_log)

# Travel distance to nearest town
RR_data$TravDist_log <- log10(RR_data$TravDist +0.1) 
RR_data$TravDist_log <- scale(RR_data$TravDist_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="TravDist") ,"Mean"]<- attr(RR_data$TravDist_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="TravDist") ,"SD"]<- attr(RR_data$TravDist_log, "scaled:scale") #sd

RR_data$TravDist_log <- as.numeric(RR_data$TravDist_log)

# Food biomass
RR_data$FoodBiomass_log <- log10(RR_data$FoodBiomass +0.1) 
RR_data$FoodBiomass_log <- scale(RR_data$FoodBiomass_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="FoodBiomass") ,"Mean"]<- attr(RR_data$FoodBiomass_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="FoodBiomass") ,"SD"]<- attr(RR_data$FoodBiomass_log, "scaled:scale") #sd
RR_data$FoodBiomass_log <- as.numeric(RR_data$FoodBiomass_log)

# Stunting prevalence
RR_data$Stunting_log <- log10(RR_data$Stunting + 0.1)
RR_data$Stunting_log <- scale(RR_data$Stunting_log, center=T, scale=T)
df_rescale[which(df_rescale$Predictor=="Stunting") ,"Mean"]<- attr(RR_data$Stunting_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="Stunting") ,"SD"]<- attr(RR_data$Stunting_log, "scaled:scale") #sd
RR_data$Stunting_log <- as.numeric(RR_data$Stunting_log)

# Human population density
RR_data$PopDens_log <- log10(RR_data$PopDens + 0.1)
RR_data$PopDens_log <- scale(RR_data$PopDens_log, center=T, scale=T)
df_rescale[which(df_rescale$Predictor=="PopDens") ,"Mean"]<- attr(RR_data$PopDens_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="PopDens") ,"SD"]<- attr(RR_data$PopDens_log, "scaled:scale") #sd
RR_data$PopDens_log <- as.numeric(RR_data$PopDens_log)

write.csv(df_rescale,"Data/Spatial/Rescaling_values.csv", row.names = F)
rm(df_rescale)


# Establish reference levels for categorical variables:
RR_data$Reserve<-as.factor(RR_data$Reserve)
RR_data$Reserve<-relevel(RR_data$Reserve, ref ="0")

RR_data$Traded<-as.factor(RR_data$Traded)
RR_data$Traded<-relevel(RR_data$Traded, ref ="No")

RR_data$Trophic_Level<-as.factor(RR_data$Trophic_Level)
RR_data$Trophic_Level<-relevel(RR_data$Trophic_Level, ref ="Omnivore")

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)

#### check multicollinearity ####
dcor <- RR_data %>% 
  dplyr::select(DistHunt_log, BodyMass_log, FoodBiomass_log, Stunting_log,
                PopDens_log, TravDist_log, HWI_log) %>% as.matrix()

corrplot::corrplot(cor(dcor), type = 'lower',addCoef.col = "black", tl.col ="#434343" ) # some high values... check vif

RR_data$log_RR <- log10(RR_data$RR +0.000001) # avoid log10(0) = -INF. Not necessary for hurdle
m <- lmer(log_RR ~ DistHunt_log + BodyMass_log + HWI_log + TravDist_log +
            PopDens_log + FoodBiomass_log + Stunting_log +(1|BirdTree_Species), data=RR_data)
performance::check_collinearity(m) # VIF< 2 for all predictors

rm(list=c("dcor", "m"))

##### Modelling hunting impacts (lognormal hurdle) ####

hurdle_formula <- bf(RR ~ DistHunt_log*BodyMass_log + I(DistHunt_log^2)*I(BodyMass_log^2) +
                        HWI_log + I(HWI_log^2) + TravDist_log + I(TravDist_log^2)+ 
                        Traded + Stunting_log + I(Stunting_log^2) + PopDens_log + 
                         I(PopDens^2) + Reserve + Trophic_Level + FoodBiomass_log + 
                       (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum),
                      hu ~ DistHunt_log*BodyMass_log + I(DistHunt_log^2)*I(BodyMass_log^2) + 
                       HWI_log + I(HWI_log^2) + TravDist_log + I(TravDist_log^2) + 
                       Traded + Stunting_log + I(Stunting_log^2) + PopDens_log + 
                       I(PopDens_log^2) + FoodBiomass_log + Reserve + Trophic_Level + 
                      (1|Dataset) + (1|CountryNum) + (1|gr(BirdTree_Species, cov = phylo_cov)))

pr <-c(prior(normal(0, 2.5), class='Intercept'), prior(normal(0, 1), class='b'))#weakly informative priors

# Import and prepare consensus tree from Jetz et al 2012. 

consensus_tree <- ape::read.nexus("Data/Phylo/ConsensusTree150_05credibility.tree")
consensus_tree$tip.label

#Create phylogenetic correlation matrix
RR_data$BirdTree_Species <- sub(" ", "_", RR_data$BirdTree_Species)
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label), ]

#Cyanoderma chrysaeum = Stachyris chrysaea
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label),"BirdTree_Species"] <- "Stachyris_chrysaea"

drops <-consensus_tree$tip.label[!consensus_tree$tip.label %in% RR_data$BirdTree_Species]
phylo_tree <- drop.tip(consensus_tree, drops)
phylo_cov <- ape::vcv.phylo(phylo_tree, corr=FALSE) #correlation matrix. If corr= TRUE => vcv matrix
rm(list=c("drops", "phylo_tree", "consensus_tree"))

tic()
hurdle_model <- brm(hurdle_formula, 
                    data = RR_data,  
                    chains=4, cores=8, iter = 4000,warmup = 2000, thin = 2, 
                    family=hurdle_lognormal(), prior=pr,  
                    data2 = list(phylo_cov = phylo_cov), seed =123,
                    control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

save(hurdle_model, file="Results/hurdle_model.RData")






=======
# Script: Bayesian modelling: Bird's ratio responses to hunting
# Author: Iago Ferreiro Arias
# Date: 3th October 2022

###############################################################################
library(brms)
library(ape)
library(spdep)
library(DHARMa)
library(ggplot2)
library(tictoc)
library(corrplot)
library(dplyr)
library(lme4)
library(performance)
################################################################################

RR_data <- read.csv("Data/Bird_RR_data.csv")

# Log10 and Scale predictors:

# Body mass
RR_data$BodyMass_log <- log10(RR_data$Body_Mass)
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)

#extract  mean and sd to rescale the same variable in the prediction dataset
# Rescaling precitors for projections: data$predictor <- (data$predictor - predictorMean)/predictorSD

df_rescale <- data.frame(Predictor=c("BodyMass", "HWI","DistHunt", "TravDist",
                                     "FoodBiomass", "Stunting", "PopDens"), 
                         Mean=NA, SD=NA)

df_rescale[which(df_rescale$Predictor=="BodyMass") ,"Mean"]<- attr(RR_data$BodyMass_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="BodyMass") ,"SD"]<- attr(RR_data$BodyMass_log, "scaled:scale") #sd

RR_data$BodyMass_log <- as.numeric(RR_data$BodyMass_log)


# Hand Wing Index
RR_data$HWI_log <- log10(RR_data$HWI)
RR_data$HWI_log <- scale(RR_data$HWI_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="HWI") ,"Mean"]<- attr(RR_data$HWI_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="HWI") ,"SD"]<- attr(RR_data$HWI_log, "scaled:scale") #sd

RR_data$HWI_log <- as.numeric(RR_data$HWI_log)


# Distance to hunter's access point

RR_data$DistHunt_log <- log10(RR_data$Dist_Hunters +0.1) 
RR_data$DistHunt_log <- scale(RR_data$DistHunt_log, center = T, scale = T)

df_rescale[which(df_rescale$Predictor=="DistHunt") ,"Mean"]<- attr(RR_data$DistHunt_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="DistHunt") ,"SD"]<- attr(RR_data$DistHunt_log, "scaled:scale") #sd

RR_data$DistHunt_log <- as.numeric(RR_data$DistHunt_log)

# Travel distance to nearest town
RR_data$TravDist_log <- log10(RR_data$TravDist +0.1) 
RR_data$TravDist_log <- scale(RR_data$TravDist_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="TravDist") ,"Mean"]<- attr(RR_data$TravDist_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="TravDist") ,"SD"]<- attr(RR_data$TravDist_log, "scaled:scale") #sd

RR_data$TravDist_log <- as.numeric(RR_data$TravDist_log)

# Food biomass
RR_data$FoodBiomass_log <- log10(RR_data$FoodBiomass +0.1) 
RR_data$FoodBiomass_log <- scale(RR_data$FoodBiomass_log, center = T, scale = T)
df_rescale[which(df_rescale$Predictor=="FoodBiomass") ,"Mean"]<- attr(RR_data$FoodBiomass_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="FoodBiomass") ,"SD"]<- attr(RR_data$FoodBiomass_log, "scaled:scale") #sd
RR_data$FoodBiomass_log <- as.numeric(RR_data$FoodBiomass_log)

# Stunting prevalence
RR_data$Stunting_log <- log10(RR_data$Stunting + 0.1)
RR_data$Stunting_log <- scale(RR_data$Stunting_log, center=T, scale=T)
df_rescale[which(df_rescale$Predictor=="Stunting") ,"Mean"]<- attr(RR_data$Stunting_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="Stunting") ,"SD"]<- attr(RR_data$Stunting_log, "scaled:scale") #sd
RR_data$Stunting_log <- as.numeric(RR_data$Stunting_log)

# Human population density
RR_data$PopDens_log <- log10(RR_data$PopDens + 0.1)
RR_data$PopDens_log <- scale(RR_data$PopDens_log, center=T, scale=T)
df_rescale[which(df_rescale$Predictor=="PopDens") ,"Mean"]<- attr(RR_data$PopDens_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="PopDens") ,"SD"]<- attr(RR_data$PopDens_log, "scaled:scale") #sd
RR_data$PopDens_log <- as.numeric(RR_data$PopDens_log)

write.csv(df_rescale,"Data/Spatial/Rescaling_values.csv", row.names = F)
rm(df_rescale)


# Establish reference levels for categorical variables:
RR_data$Reserve<-as.factor(RR_data$Reserve)
RR_data$Reserve<-relevel(RR_data$Reserve, ref ="0")

RR_data$Traded<-as.factor(RR_data$Traded)
RR_data$Traded<-relevel(RR_data$Traded, ref ="No")

RR_data$Trophic_Level<-as.factor(RR_data$Trophic_Level)
RR_data$Trophic_Level<-relevel(RR_data$Trophic_Level, ref ="Omnivore")

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)

#### check multicollinearity ####
dcor <- RR_data %>% 
  dplyr::select(DistHunt_log, BodyMass_log, FoodBiomass_log, Stunting_log,
                PopDens_log, TravDist_log, HWI_log) %>% as.matrix()

corrplot::corrplot(cor(dcor), type = 'lower',addCoef.col = "black", tl.col ="#434343" ) # some high values... check vif

RR_data$log_RR <- log10(RR_data$RR +0.000001) # avoid log10(0) = -INF. Not necessary for hurdle
m <- lmer(log_RR ~ DistHunt_log + BodyMass_log + HWI_log + TravDist_log +
            PopDens_log + FoodBiomass_log + Stunting_log +(1|BirdTree_Species), data=RR_data)
performance::check_collinearity(m) # VIF< 2 for all predictors

rm(list=c("dcor", "m"))

##### Modelling hunting impacts (lognormal hurdle) ####

hurdle_formula <- bf(RR ~ DistHunt_log*BodyMass_log + I(DistHunt_log^2)*I(BodyMass_log^2) +
                        HWI_log + I(HWI_log^2) + TravDist_log + I(TravDist_log^2)+ 
                        Traded + Stunting_log + I(Stunting_log^2) + PopDens_log + 
                         I(PopDens^2) + Reserve + Trophic_Level + FoodBiomass_log + 
                       (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum),
                      hu ~ DistHunt_log*BodyMass_log + I(DistHunt_log^2)*I(BodyMass_log^2) + 
                       HWI_log + I(HWI_log^2) + TravDist_log + I(TravDist_log^2) + 
                       Traded + Stunting_log + I(Stunting_log^2) + PopDens_log + 
                       I(PopDens_log^2) + FoodBiomass_log + Reserve + Trophic_Level + 
                      (1|Dataset) + (1|CountryNum) + (1|gr(BirdTree_Species, cov = phylo_cov)))

pr <-c(prior(normal(0, 2.5), class='Intercept'), prior(normal(0, 1), class='b'))#weakly informative priors

# Import and prepare consensus tree from Jetz et al 2012. 

consensus_tree <- ape::read.nexus("Data/Phylo/ConsensusTree150_05credibility.tree")
consensus_tree$tip.label

#Create phylogenetic correlation matrix
RR_data$BirdTree_Species <- sub(" ", "_", RR_data$BirdTree_Species)
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label), ]

#Cyanoderma chrysaeum = Stachyris chrysaea
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label),"BirdTree_Species"] <- "Stachyris_chrysaea"

drops <-consensus_tree$tip.label[!consensus_tree$tip.label %in% RR_data$BirdTree_Species]
phylo_tree <- drop.tip(consensus_tree, drops)
phylo_cov <- ape::vcv.phylo(phylo_tree, corr=FALSE) #correlation matrix. If corr= TRUE => vcv matrix
rm(list=c("drops", "phylo_tree", "consensus_tree"))

tic()
hurdle_model <- brm(hurdle_formula, 
                    data = RR_data,  
                    chains=4, cores=8, iter = 4000,warmup = 2000, thin = 2, 
                    family=hurdle_lognormal(), prior=pr,  
                    data2 = list(phylo_cov = phylo_cov), seed =123,
                    control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

save(hurdle_model, file="Results/hurdle_model.RData")






>>>>>>> Stashed changes
