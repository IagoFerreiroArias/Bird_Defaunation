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
library(ggpubr)
################################################################################

RR_data <- read.csv("Data/Bird_RR_data.csv")

length(unique(RR_data$Study))
length(unique(RR_data$Species))

RR_data$RR_log <- log(RR_data$RR)
summary(RR_data[which(RR_data$RR_log !='-Inf'),]$RR_log)
summary(RR_data$RR)
hist(RR_data$RR_log)

#percentage of local extirpations
length(which(RR_data$bin == 0)) / length(RR_data$bin)* 100

summary(RR_data$RR)

length(RR_data[which(RR_data$RR>5),"RR"]) #21 observations

# Categorical predictor: non-hunted species vs species hunted (for food and/or trade)
RR_data$Hunted <- ifelse(RR_data$Food == "No" & RR_data$Traded == "No", "No", "Yes")
unique(RR_data[RR_data$Hunted=="No","BirdTree_Species"])

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
table(RR_data$Realm)

RR_data$Bin <- ifelse(RR_data$RR==0, "Local extirpation",
                      ifelse(RR_data$RR==1, "No changes in abundance", 
                             ifelse(RR_data$RR>1, "Increase in abundance", "Declines in abundance")))

RR_data$Bin <- as.factor(RR_data$Bin)
RR_data$Bin <-factor(RR_data$Bin, levels = c("Local extirpation", "Declines in abundance","No changes in abundance", "Increase in abundance"))
levels(RR_data$Bin)

d1 <- ggplot(data=RR_data, aes(RR, colour=Bin, fill=Bin)) + geom_histogram(bins= 100) + 
  theme_bw() + ylab("Count") + theme(legend.position = c(0.8,0.9), legend.background = element_blank()) +
  scale_color_manual(values=c("#D40404", "#FAA37A","#FAE77A", "#98BF99"))+ geom_vline(xintercept=1, linetype="dashed") +
  labs(fill='Response to hunting pressure', colour= "Response to hunting pressure") + 
  scale_fill_manual(values=c("#D40404", "#FAA37A","#FAE77A", "#98BF99")) + xlim(-0.3, max(RR_data$RR)) 

RR_data$Bin <- ifelse(RR_data$RR_log==0, "No changes in abundance",
                      ifelse(RR_data$RR>1, "Increase in abundance",  "Reduction in abundance"))

RR_data$Bin <- as.factor(RR_data$Bin)
RR_data$Bin <-factor(RR_data$Bin, levels = c("Reduction of abundance","No changes in abundance", "Increase in abundance"))
d2<- ggplot(data=RR_data[which(RR_data$RR_log !='-Inf'),], aes(RR_log, colour=Bin, fill=Bin)) + geom_histogram(bins= 500) + 
  theme_bw() + ylab("Count") + theme(legend.position = c(0.8,0.8), legend.background = element_blank()) +
  scale_color_manual(values=c("#FAA37A","#FAE77A","#98BF99"))+ geom_vline(xintercept=0, linetype="dashed") +
  labs(fill='Response to hunting pressure', colour= "Response to hunting pressure") + 
  scale_fill_manual(values=c("#FAA37A","#FAE77A","#98BF99"))  + xlab("log(RR)") + xlim(-5,5)

ggpubr::ggarrange(d1,d2)

rm(list=c("d1","d2"))

RR_data$Size <- ifelse(RR_data$Body_Mass < 30, "Small",
                       ifelse(RR_data$Body_Mass > 1500, "Large", "Medium"))

RR_data$Dist <- ifelse(RR_data$Dist_Hunters  < 3, "Short",
                       ifelse(RR_data$Dist_Hunters  > 8, "Large", "Medium"))

RR_data$bin2 <- ifelse(RR_data$bin==0, "Extirpated", "Non extirpated")
RR_data$reserve2 <- ifelse(RR_data$Reserve==0, "Outside", "Inside")
table(RR_data$bin2, RR_data$reserve2)

#extirpations <- subset(RR_data, bin==0)
#non_extirpations <- subset(RR_data, bin==1)
#outside <- subset(non_extirpations, Reserve==0)
#inside <- subset(non_extirpations, Reserve==1)
#wgs84<-CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs")

#outside <- SpatialPointsDataFrame(coords=cbind(outside$Longitude,outside$Latitude),
#                                       data=outside,proj4string = wgs84)
#inside <- SpatialPointsDataFrame(coords=cbind(inside$Longitude,inside$Latitude),
#                                  data=inside,proj4string = wgs84)

#Reserve <- raster("Data/Spatial/Projections/Reserves.tif")

#mapview::mapview(Reserve, col="darkgreen", maxpixels =  5000000) + 
#  mapview::mapview(outside, col.regions="darkred") + 
#  mapview::mapview(inside, col.regions="darkblue")

# Log10 and Scale predictors:

# Body mass
RR_data$BodyMass_log <- log10(RR_data$Body_Mass/1000)
RR_data$BodyMass_log <- scale(RR_data$BodyMass_log, center = T, scale = T)

#extract  mean and sd to rescale the same variable in the prediction dataset
# Rescaling precitors for projections: data$predictor <- (data$predictor - predictorMean)/predictorSD

df_rescale <- data.frame(Predictor=c("BodyMass", "DistHunt", "TravDist",
                                     "Stunting", "PopDens", "NPP"), 
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

# NPP
RR_data$NPP_log <- log10(RR_data$NPP + 0.1)
RR_data$NPP_log <- scale(RR_data$NPP_log, center=T, scale=T)
df_rescale[which(df_rescale$Predictor=="NPP") ,"Mean"]<- attr(RR_data$NPP_log, "scaled:center") #mean
df_rescale[which(df_rescale$Predictor=="NPP") ,"SD"]<- attr(RR_data$NPP_log, "scaled:scale") #sd
RR_data$NPP_log <- as.numeric(RR_data$NPP_log)

write.csv(df_rescale,"Results/Rescaling_values.csv", row.names = F)
#rm(df_rescale)

# Establish reference levels for categorical variables:
RR_data$Reserve<-as.factor(RR_data$Reserve)
RR_data$Reserve<-relevel(RR_data$Reserve, ref ="0")

RR_data$BirdTree_Species <- as.factor(RR_data$BirdTree_Species)
RR_data$Dataset <- as.factor(RR_data$Dataset)
RR_data$CountryNum <- as.factor(RR_data$CountryNum)

summary(RR_data$NPP_log)
RR_data$NPP_cat <- ifelse(RR_data$NPP_log>0.77, "High",
                          ifelse(RR_data$NPP_log> -0.7,"Low", "Moderate"))
#ggplot(RR_data, aes(x=Dist_Hunters, fill=NPP_cat, alpha=0.5)) + geom_density() + theme_bw() + facet_grid(.~Size)
#ggplot(RR_data, aes(x=Body_Mass, fill=NPP_cat, alpha=0.5)) + geom_density() + theme_bw()+ xlim(0,5000)
#RR_data %>% group_by(NPP_cat) %>% count()

RR_data <- subset(RR_data, RR<5) # remove outliers

#### check multicollinearity ####
dcor <- RR_data %>% 
  dplyr::select(DistHunt_log, BodyMass_log, Stunting_log, 
                PopDens_log, TravDist_log, NPP_log) %>% as.matrix()

corrplot::corrplot(cor(dcor), type = 'lower',addCoef.col = "black", tl.col ="#434343" ) # some high values... check vif

RR_data$log_RR <- log(RR_data$RR +0.000001) # avoid log10(0) = -INF. Not necessary for hurdle

m <- lmer(log_RR ~ DistHunt_log + BodyMass_log + TravDist_log + 
            PopDens_log + Stunting_log + NPP_log + (1|BirdTree_Species), data=RR_data)
vif <- plot(performance::check_collinearity(m)) # VIF< 2 for all predictors

rm(list=c("dcor", "m", "vif"))

##### Modelling hunting impacts (lognormal hurdle) ####

pr <-c(prior(normal(0, 10), class='Intercept'), prior(normal(0, 0.5), class='b'))#weakly informative priors

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
rm(list=c("drops", "phylo_tree", "consensus_tree"))

# hurdle generalized mixed model
hurdle_formula <- bf(RR ~ DistHunt_log + BodyMass_log + Stunting_log + TravDist_log + 
                       PopDens_log +  NPP_log + Reserve + 
                       BodyMass_log:DistHunt_log + BodyMass_log:TravDist_log + 
                       BodyMass_log:Stunting_log + NPP_log:DistHunt_log  +
                       (1|gr(BirdTree_Species, cov = phylo_cov)) + (1|CountryNum),
                     hu ~ DistHunt_log + BodyMass_log + Stunting_log + TravDist_log + 
                       PopDens_log +  NPP_log + Reserve + 
                       BodyMass_log:DistHunt_log + BodyMass_log:TravDist_log + 
                       BodyMass_log:Stunting_log + NPP_log:DistHunt_log +
                       (1|CountryNum) + (1|gr(BirdTree_Species, cov = phylo_cov)))


tic()
hurdle_model <- brm(hurdle_formula, 
                  data = RR_data,  
                  chains=4, cores=4, iter = 4000,warmup = 2000, thin = 2, 
                  family=hurdle_lognormal(), prior=pr,  
                  data2 = list(phylo_cov = phylo_cov), seed=123,
                  control = list(adapt_delta =0.99, max_treedepth=15))
toc() #11817.4 sec elapsed

save(hurdle_model, file="Results/hurdle_model.RData")
