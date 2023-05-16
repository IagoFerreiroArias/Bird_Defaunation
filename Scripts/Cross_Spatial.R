# Script: Spatial cross-validation for hurdle model
# Author: Iago Ferreiro Arias
# Date: 12th May 2023

library(brms) #bayessian hurdle modelling 
library(loo) # k-fold crossvalidation
library(ape) # vcov matrix for hurdle models
library(tictoc)
library(raster) #create spatial blocks
library(DHARMa) #Spatial autocorrelation


########### Recreate data (cov matrix and data) for model (same code as RR_BRMS.R) ###########

RR_data <- read.csv("Data/Bird_RR_data.csv") # import clean dataset
RR_data$RR_log <- log(RR_data$RR) # check log transformation

RR_data$BodyMass_log <- log10(RR_data$Body_Mass/1000)#log10 transformation and into kg
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)#scale and center predictor
RR_data$BodyMass_log <- as.numeric(RR_data$BodyMass_log)#turn into numeric vector

RR_data$DistHunt_log <- log10(RR_data$Dist_Hunters +0.1) #log10 transformation
RR_data$DistHunt_log <- scale(RR_data$DistHunt_log, center = T, scale = T)#scale and center predictor
RR_data$DistHunt_log <- as.numeric(RR_data$DistHunt_log)#turn into numeric vector

RR_data$TravDist_log <- log10(RR_data$TravDist +0.1) #log10 transformation
RR_data$TravDist_log <- scale(RR_data$TravDist_log, center = T, scale = T)
RR_data$TravDist_log <- as.numeric(RR_data$TravDist_log)#turn into numeric vector

RR_data$Stunting_log <- log10(RR_data$Stunting + 0.1)#log10 transformation
RR_data$Stunting_log <- scale(RR_data$Stunting_log, center=T, scale=T)#scale and center predictor
RR_data$Stunting_log <- as.numeric(RR_data$Stunting_log)#turn into numeric vector


RR_data$PopDens_log <- log10(RR_data$PopDens + 0.1) #log10 transformation
RR_data$PopDens_log <- scale(RR_data$PopDens_log, center=T, scale=T) #scale and center predictor
RR_data$PopDens_log <- as.numeric(RR_data$PopDens_log) #turn into numeric vector


RR_data$Reserve<-as.factor(RR_data$Reserve) 
RR_data$Reserve<-relevel(RR_data$Reserve, ref ="0")# Establish reference levels for categorical variables

RR_data$Traded<-as.factor(RR_data$Traded)
RR_data$Traded<-relevel(RR_data$Traded, ref ="No")# Establish reference levels for categorical variables

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)



consensus_tree <- ape::read.nexus("Data/Phylo/ConsensusTree150_05credibility.tree") # Import and prepare consensus tree from Jetz et al 2012. 
consensus_tree$tip.label
RR_data$BirdTree_Species <- sub(" ", "_", RR_data$BirdTree_Species)
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label), ]
RR_data[-which(RR_data$BirdTree_Species %in% consensus_tree$tip.label),"BirdTree_Species"] <- "Stachyris_chrysaea" #Cyanoderma chrysaeum = Stachyris chrysaea
drops <-consensus_tree$tip.label[!consensus_tree$tip.label %in% RR_data$BirdTree_Species]
phylo_tree <- drop.tip(consensus_tree, drops)
phylo_cov <- ape::vcv.phylo(phylo_tree, corr=FALSE) # correlation matrix. If corr= TRUE => vcv matrix
rm(list=c("drops", "phylo_tree", "consensus_tree"))



pr <-c(prior(normal(0, 10), class='Intercept'), prior(normal(0, 1), class='b'))#weakly informative priors
hurdle_formula <- bf(RR ~ s(DistHunt_log,k=3) + s(BodyMass_log, k=3) +  s(TravDist_log, k=3)+
                       s(Stunting_log, k=3) + s(PopDens_log, k=3) + Traded + Reserve +
                       (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|Dataset) + (1|CountryNum),
                     hu ~  s(DistHunt_log, k=3) + s(BodyMass_log, k=3) +  s(TravDist_log, k=3)+
                       s(Stunting_log, k=3) + s(PopDens_log, k=3) + Reserve + Traded +
                       (1|Dataset) + (1|CountryNum) + (1|gr(BirdTree_Species, cov = phylo_cov)))


# Spatial-block cross-validation for hurdle model

blockSize <- 5 #degrees
Xvec<-seq(-180, 180, by=blockSize)
Yvec<-seq(-38, 43, by=blockSize) #Degrees of latitude forest extension: - 37.44 , 42.60 (gis layer)
grid<-expand.grid(Longitude=Xvec, Latitude=Yvec)
grid$ID<-1:nrow(grid)
SpatialBlocks<-raster::rasterFromXYZ(grid)
RR_data$SpatialBlock<-raster::extract(SpatialBlocks, RR_data[,c('Longitude','Latitude')])

rmse<-rep(NA, 10)

RR_data$foldsSpatial <- kfold_split_grouped(K=10, x=RR_data$SpatialBlock)

for (i in unique(RR_data$foldsSpatial)){
  
  testing<-RR_data[RR_data$foldsSpatial %in% i,]
  training<-RR_data[!RR_data$foldsSpatial %in% i,]
  Traded<-all(unique(testing$Traded) %in% unique(training$Traded)) #always true?
  Reserve<-all(unique(testing$Reserve) %in% unique(training$Reserve)) #always true?
  if(!all(Traded, Reserve)) {next}
  
  mod <- brm(hurdle_formula, 
             data = RR_data,  
             chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
             family=hurdle_lognormal(), prior=pr,  
             data2 = list(phylo_cov = phylo_cov), seed=123,
             control = list(adapt_delta =0.99, max_treedepth=15))
  
  pred<-predict(mod, newdata=testing, summary=TRUE, re_formula=NA)
  rmse[i]<-sqrt(mean((pred[,1] - testing$RR)^2))
  
}


rmseSpatial<-rmse
write.csv(rmseSpatial, 'Results/Spatial_CrossValidation.csv', row.names=FALSE)


# Check spatial autocorrelation model residuals

RR_data$XY<-paste0(RR_data$Longitude,'_',RR_data$Latitude)
XY<-unique(RR_data$XY)
XY2<-matrix(as.numeric(unlist(strsplit(XY, '_'))), ncol=2, byrow=TRUE)
X<-XY2[,1]
Y<-XY2[,2]

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(hurdle_gam)),
  observedResponse = RR_data$RR,
  fittedPredictedResponse = apply(t(posterior_epred(hurdle_gam)), 1, mean),
  integerResponse = TRUE)

modelcheck_SP<-recalculateResiduals(model_check, group=RR_data$XY)
MoranI<-testSpatialAutocorrelation(modelcheck_SP, x=X, y=Y) #p>0.05 means no spatial autocorrelation
MoranI<-data.frame(Observed=MoranI$statistic[1], Expected=MoranI$statistic[2], SD=MoranI$statistic[3], Pvalue=MoranI$p.value)
row.names(MoranI)<- "MoranI"
write.csv(MoranI,'Results/MoranI.csv', row.names=FALSE)


rm(list=c("blockSize", "i", "rmse", "X", "Xvec", "Y", "Yvec", "XY", "XY2",
          "MoranI", "modelcheck_SP", "model_check", "SpatialBlocks", "grid",
          "Reserve", "Traded"))
