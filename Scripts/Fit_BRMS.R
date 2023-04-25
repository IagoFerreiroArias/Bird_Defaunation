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
?log
RR_data$RR_log <- log(RR_data$RR)
summary(RR_data[which(RR_data$RR_log !='-Inf'),"RR_log"])
hist(RR_data[which(RR_data$RR_log !='-Inf'),"RR_log"], xlab="RR", main ="RR")


# Log10 and Scale predictors:

# Body mass
RR_data$BodyMass_log <- log10(RR_data$Body_Mass)
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)

#extract  mean and sd to rescale the same variable in the prediction dataset
# Rescaling precitors for projections: data$predictor <- (data$predictor - predictorMean)/predictorSD

df_rescale <- data.frame(Predictor=c("BodyMass", "DistHunt", "TravDist",
                                     "Stunting", "PopDens"), 
                         Mean=NA, SD=NA)

df_rescale[which(df_rescale$Predictor=="BodyMass") ,"Mean"]<- attr(RR_data$BodyMass_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="BodyMass") ,"SD"]<- attr(RR_data$BodyMass_log, "scaled:scale") #sd

RR_data$BodyMass_log <- as.numeric(RR_data$BodyMass_log)


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

# Prevalence stunting 
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
RR_data$Trophic_Level<-relevel(RR_data$Trophic_Level, ref ="Herbivore")

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)

#### check multicollinearity ####
dcor <- RR_data %>% 
  dplyr::select(DistHunt_log, BodyMass_log,  Stunting_log,
                PopDens_log, TravDist_log) %>% as.matrix()

corrplot::corrplot(cor(dcor), type = 'lower',addCoef.col = "black", tl.col ="#434343" ) # some high values... check vif

RR_data$log_RR <- log10(RR_data$RR +0.000001) # avoid log10(0) = -INF. Not necessary for hurdle

m <- lmer(log_RR ~ DistHunt_log + BodyMass_log + TravDist_log +
            PopDens_log + Stunting_log +(1|BirdTree_Species), data=RR_data)
performance::check_collinearity(m) # VIF< 2 for all predictors

rm(list=c("dcor", "m"))

# Check distribution of body mass in function of trade category and trophic level

#ggpubr::ggarrange(ggplot(RR_data) + geom_density(aes(log10(Body_Mass/1000), fill = Traded, alpha=0.8)) +
#                    ylab ("Density") + xlab("log10(Body Mass) (kg)") + theme_bw() + scale_alpha(guide = 'none')+
#                    scale_fill_manual(values=c("#DAF7A6", "#FF5733"),"Traded")+ 
#                    theme(legend.position=c(0.9,0.5), legend.background = element_blank()),
#                  ggplot(RR_data) + geom_density(aes(log10(Body_Mass/1000), fill = Trophic_Level, alpha=0.9)) +
#                    ylab ("Density") + xlab("log10(Body Mass) (kg)") + theme_bw() + scale_alpha(guide = 'none')+
#                    theme(legend.position=c(0.9,0.5), legend.background = element_blank())+ 
#                    scale_fill_manual(values=c("#DAF7A6", "#FF5733", "#FFC300"),"Trophic level")) 


##### Modelling hunting impacts (changes in abundance and probability of extirpation) ####


gaussian_formula <- bf(RR_log ~ s(DistHunt_log, BodyMass_log, k=3) + s(TravDist_log, k=3) +
                         s(Stunting_log, k=3) + s(PopDens_log, k=3) + Reserve + Traded + 
                         (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))


bin_formula <- bf(bin ~ s(DistHunt_log, BodyMass_log, k=3) + s(TravDist_log, k=3) + 
                    s(Stunting_log, k=3) + s(PopDens_log , k=3) + Traded + Reserve + 
                    (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))

pr <-c(prior(normal(0, 1), class='Intercept'), prior(normal(0, 1), class='b'))#weakly informative priors

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
gaussian <- brm(gaussian_formula, 
                    data = RR_data[which(RR_data$RR_log !='-Inf'),],  
                    chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                    family="gaussian", prior=pr,  
                    data2 = list(phylo_cov = phylo_cov), seed=123,
                    control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

save(gaussian, file="Results/gaussian_model.RData")


tic()
binomial <- brm(bin_formula, 
                data = RR_data,  
                chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                family="binomial", prior=pr,  
                data2 = list(phylo_cov = phylo_cov), seed=123,
                control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

formula(binomial)

save(binomial, file="Results/binomial_model.RData")

##

load("Results/binomial_model.Rdata")
load("Results/gaussian_model.Rdata")

summary(binomial)
summary(gaussian)

r2_bayes(binomial)
r2_bayes(gaussian)
pp_check(binomial, ndraws = 100) + xlim(0,2)
pp_check(gaussian, ndraws = 100) + xlim(0,3)

#### DHARMa diagnosis ####

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(binomial)),
  observedResponse = RR_data$bin,
  fittedPredictedResponse = apply(t(posterior_epred(binomial)), 1, mean),
  integerResponse = TRUE)


model_check2 <- createDHARMa(
  simulatedResponse = t(posterior_predict(gaussian)),
  observedResponse = RR_data[which(RR_data$RR_log !='-Inf'),"RR_log"],
  fittedPredictedResponse = apply(t(posterior_epred(gaussian)), 1, mean),
  integerResponse = TRUE)

plot(model_check2)


##### Modelling hunting impacts (changes in abundance and probability of extirpation) ####

# Body mass x dist sin traded

gaussian_formula1 <- bf(RR_log ~ s(DistHunt_log, BodyMass_log, k=3) + s(TravDist_log, k=3) +
                         s(Stunting_log, k=3) + s(PopDens_log , k=3) + Reserve + 
                         (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))


bin_formula1 <- bf(bin ~ s(DistHunt_log, BodyMass_log, k=3) + s(TravDist_log, k=3) + 
                    s(Stunting_log, k=3) + s(PopDens_log , k=3) + Reserve + 
                    (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))

tic()
gaussian1 <- brm(gaussian_formula1, 
                data = RR_data[which(RR_data$RR_log !='-Inf'),],  
                chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                family="gaussian", prior=pr,  
                data2 = list(phylo_cov = phylo_cov), seed=123,
                control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

summary(gaussian1)

load("Results/gaussian_model1.RData")
conditional_effects(gaussian1, method="predict", prob=0.95)#

save(gaussian1, file="Results/gaussian_model1.RData")

tic()
binomial1 <- brm(bin_formula1, 
                data = RR_data,  
                chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                family="binomial", prior=pr,  
                data2 = list(phylo_cov = phylo_cov), seed=123,
                control = list(adapt_delta =0.95, max_treedepth=15))

load("Results/binomial_model1.RData")
conditional_effects(binomial1, method="predict", prob=0.95)#
bayes_R2(binomial1)
toc() #11817.4 sec elapsed

summary(binomial1)


save(binomial1, file="Results/binomial_model1.RData")

conditional_effects(binomial1, method="predict", prob=0.95)#

# traded x dist sin traded

gaussian_formula2 <- bf(RR_log ~ s(DistHunt_log, by=Traded, k=3) + s(TravDist_log, k=3) +
                          s(Stunting_log, k=3) + s(PopDens_log , k=3) + Reserve + 
                          (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))


bin_formula2 <- bf(bin ~ s(DistHunt_log,by=Traded, k=3) + s(TravDist_log, k=3) + 
                     s(Stunting_log, k=3) + s(PopDens_log , k=3) + Reserve + 
                     (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))

tic()
gaussian2 <- brm(gaussian_formula2, 
                 data = RR_data[which(RR_data$RR_log !='-Inf'),],  
                 chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                 family="gaussian", prior=pr,  
                 data2 = list(phylo_cov = phylo_cov), seed=123,
                 control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed


save(gaussian2, file="Results/gaussian_model2.RData")
load("Results/gaussian_model2.RData")

conditional_effects(gaussian2, method="predict", prob=0.95)#

tic()
binomial2 <- brm(bin_formula2, 
                 data = RR_data,  
                 chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                 family="binomial", prior=pr,  
                 data2 = list(phylo_cov = phylo_cov), seed=123,
                 control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed


save(binomial2, file="Results/binomial_model2.RData")
load("Results/binomial_model2.RData")

conditional_effects(binomial2, method="predict", prob=0.95)#

#  sin travel time

gaussian_formula3 <- bf(RR_log ~ s(DistHunt_log, BodyMass_log, k=3) + 
                          s(Stunting_log, k=3) + s(PopDens_log, k=3) + Reserve + 
                          (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))


bin_formula3 <- bf(bin ~ s(DistHunt_log,BodyMass_log, k=3) + 
                     s(Stunting_log, k=3) + s(PopDens_log , k=3) + Reserve + 
                     (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum))

tic()
gaussian3 <- brm(gaussian_formula3, 
                 data = RR_data[which(RR_data$RR_log !='-Inf'),],  
                 chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                 family="gaussian", prior=pr,  
                 data2 = list(phylo_cov = phylo_cov), seed=123,
                 control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed

save(gaussian3, file="Results/gaussian_model3.RData")
load("Results/gaussian_model3.RData")

conditional_effects(gaussian3, method="predict", prob=0.95)#


tic()
binomial3 <- brm(bin_formula3, 
                 data = RR_data,  
                 chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                 family="binomial", prior=pr,  
                 data2 = list(phylo_cov = phylo_cov), seed=123,
                 control = list(adapt_delta =0.95, max_treedepth=15))

toc() #11817.4 sec elapsed


save(binomial3, file="Results/binomial_model3.RData")



