#Script: Model diagnosis, plots and cross-validations
#Author: Iago Ferreiro Arias
#Date: 16th March, 2023

library(brms)
library(marginaleffects)
library(sjPlot)
library(DHARMa)
library(performance)
library(ggplot2)
library(spdep)
library(emmeans)
library(tidybayes)

?conditional_effects
# import hurdle model

load("Results/hurdle_model.RData")

hist(log10(RR_data$Body_Mass/1000))
summary(RR_data$Body_Mass/1000)
summary(hurdle_model)


#### Posterior predictive checks ####
formula(hurdle_model)
pp_check(hurdle_model, ndraws = 100) + xlim(0,30)

#### DHARMa diagnosis ####

model_check <- createDHARMa(
  simulatedResponse = t(posterior_predict(hurdle_model)),
  observedResponse = RR_data$RR,
  fittedPredictedResponse = apply(t(posterior_epred(hurdle_model)), 1, mean),
  integerResponse = TRUE)

plot(model_check)
rm(model_check)

summary(hurdle_model)
##### Plotting effects #####

c_eff <- conditional_effects(hurdle_model, method="predict", dpar="mu", prob=0.95)#non-zero part
c_eff2 <- conditional_effects(hurdle_model, method="predict", dpar="mu", prob=0.8)
c_eff3 <- conditional_effects(hurdle_model, method="predict", dpar="mu", prob=0.5)

mu1 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[1]] +
  theme_bw()+ ylab("RR") + xlab("Distance to hunter's access point (log)") +
  geom_ribbon(data = c_eff2$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")

mu2 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[2]] +
  theme_bw()+ ylab("RR") + xlab("Body mass (log)") +
  geom_ribbon(data = c_eff2$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(mu2)

mu3 <- plot(c_eff, plot=FALSE, list(alpha=0.2))[[9]] +theme_bw() + ylab("RR") + 
  xlab("Distance to hunter's access point (log)") +  ylim(0,2) + 
  scale_fill_manual(values=c("#A10000", "#D0A277", "#FDAF17"),"Body mass")+
  scale_colour_manual(values=c("#820303", "#BC926B", "#DA9611"), "Body mass") + 
  theme(legend.position=c(0.5,0.9), legend.direction = "horizontal")

plot(mu3)

mu4 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[3]] +
  theme_bw()+ ylab("RR") + xlab("Travel time to major cities (log)") +
  geom_ribbon(data = c_eff2$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(mu4)

mu5 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[5]] +
  theme_bw()+ ylab("RR") + xlab("Prevalence of stunting (log)") +
  geom_ribbon(data = c_eff2$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(mu5)

mu6 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[6]] +
  theme_bw()+ ylab("RR") + xlab("Human population density (log)") +
  geom_ribbon(data = c_eff2$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(mu6)

mu7 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[8]] +
  theme_bw()+ ylab("RR") + xlab("Food biomass (log)") +
  geom_ribbon(data = c_eff2$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(mu7)

mu8 <- plot(c_eff, plot=FALSE)[[7]] + theme_bw() + ylab("RR")
mu9 <- plot(c_eff, plot=FALSE)[[4]] + theme_bw() + ylab("RR")



ggpubr::ggarrange(mu1,mu2, mu3, mu4, mu5, mu6, mu7, mu8 , mu9) 
ggsave(plot1, filename="Figures/Model/RR_Marginal_effects.pdf", width=400, height=313, units="mm", device="pdf")
rm(list=c("mu1", "mu2", "mu3","mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "plot1"))

### plot probability of being locally extirpated

c_eff <- conditional_effects(hurdle_model, method="fitted", dpar="hu", prob=0.95)#zero part
c_eff2 <- conditional_effects(hurdle_model, method="fitted", dpar="hu", prob=0.8)
c_eff3 <- conditional_effects(hurdle_model, method="fitted", dpar="hu", prob=0.5)

hu1 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[1]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Distance to hunter's access point (log)") +
  geom_ribbon(data = c_eff2$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu1)

hu2 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[2]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Body mass (log)") +
  geom_ribbon(data = c_eff2$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu2)

hu3 <- plot(c_eff, plot=FALSE, list(alpha=0.2))[[9]] +theme_bw() + ylab("Probability of being locally extirpated") + 
  xlab("Distance to hunter's access point (log)") + ymin(0,2) +
  scale_fill_manual(values=c("#A10000", "#D0A277", "#FDAF17"),"Body mass")+
  scale_colour_manual(values=c("#820303", "#BC926B", "#DA9611"), "Body mass") + 
  theme(legend.position=c(0.5,0.9), legend.direction = "horizontal")
plot(hu3)

hu4 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[3]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Travel time to major cities (log)") +
  geom_ribbon(data = c_eff2$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu4)

hu5 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[5]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Prevalence of stunting (log)") +
  geom_ribbon(data = c_eff2$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu5)

hu6 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[6]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Human population density (log)") +
  geom_ribbon(data = c_eff2$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu6)

hu7 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[8]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Food biomass (log)") +
  geom_ribbon(data = c_eff2$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(hu7)

hu8 <- plot(c_eff, plot=FALSE)[[7]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(hu8)
hu9 <- plot(c_eff, plot=FALSE)[[4]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(hu9)
hu10 <- plot(c_eff, plot=FALSE)[[10]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(hu10)

plot2 <- ggpubr::ggarrange(hu1,hu2, hu3, hu4, hu5, hu6, hu7, hu8 , hu9, hu10) 
ggsave(plot2, filename="Figures/Model/HU_Marginal_effects.pdf", width=400, height=313, units="mm", device="pdf")
rm(list=c("hu1", "hu2", "hu3","hu4", "hu5", "hu6", "hu7", "hu8", "hu9", "hu10"))
rm(plot2)


summary(hurdle_model)
brms::bayes_R2(hurdle_model)
#probability of being locally extirpated due to overhunting
conditional_effects(hurdle_model,dpar="hu")

formula(hurdle_model)
#data_term <- conditional_effects(model)$term
#data_term <- conditional_effects(model)$`term1:term2`


##### Plotting effects: binomial model #####
c_eff <- conditional_effects(binomial, method="fitted", prob=0.95)#non-zero part
c_eff2 <- conditional_effects(binomial, method="fitted", prob=0.8)
c_eff3 <- conditional_effects(binomial, method="fitted", prob=0.5)

# reference level in binomial response = 1, change to 0 to better 
# interpret the probability of being locally extirpated
c_eff$DistHunt_log$estimate__ <- 1 - c_eff$DistHunt_log$estimate__
c_eff$DistHunt_log$lower__ <- 1 - c_eff$DistHunt_log$lower__ 
c_eff$DistHunt_log$upper__ <- 1 - c_eff$DistHunt_log$upper__

c_eff2$DistHunt_log$lower__ <- 1 - c_eff2$DistHunt_log$lower__ 
c_eff2$DistHunt_log$upper__ <- 1 - c_eff2$DistHunt_log$upper__

c_eff3$DistHunt_log$lower__ <- 1 - c_eff3$DistHunt_log$lower__ 
c_eff3$DistHunt_log$upper__ <- 1 - c_eff3$DistHunt_log$upper__


bin1 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[1]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Distance to hunter's access point (log)") +
  geom_ribbon(data = c_eff2$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")

######### BODY MASS ######
c_eff$BodyMass_log$estimate__ <- 1 - c_eff$BodyMass_log$estimate__
c_eff$BodyMass_log$lower__ <- 1 - c_eff$BodyMass_log$lower__ 
c_eff$BodyMass_log$upper__ <- 1 - c_eff$BodyMass_log$upper__

c_eff2$BodyMass_log$lower__ <- 1 - c_eff2$BodyMass_log$lower__ 
c_eff2$BodyMass_log$upper__ <- 1 - c_eff2$BodyMass_log$upper__

c_eff3$BodyMass_log$lower__ <- 1 - c_eff3$BodyMass_log$lower__ 
c_eff3$BodyMass_log$upper__ <- 1 - c_eff3$BodyMass_log$upper__

bin2 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[2]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Body mass (log)") +
  geom_ribbon(data = c_eff2$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(bin2)

### Interaction term: distance + body mass ####
c_eff$`DistHunt_log:BodyMass_log`$estimate__ <- 1 - c_eff$`DistHunt_log:BodyMass_log`$estimate__
c_eff$`DistHunt_log:BodyMass_log`$lower__ <- 1 - c_eff$`DistHunt_log:BodyMass_log`$lower__ 
c_eff$`DistHunt_log:BodyMass_log`$upper__ <- 1 -  c_eff$`DistHunt_log:BodyMass_log`$upper__ 

bin10 <- plot(c_eff, plot=FALSE, list(alpha=0.2))[[10]] +theme_bw() + 
  ylab("Probability of being locally extirpated") + 
  xlab("Distance to hunter's access point (log)") + 
  scale_fill_manual(values=c("#A10000", "#D0A277", "#FDAF17"),"Body mass")+
  scale_colour_manual(values=c("#820303", "#BC926B", "#DA9611"), "Body mass") + 
  theme(legend.position=c(0.5,0.95), legend.direction = "horizontal")

plot(bin10)

######### Travel time ######

c_eff$TravDist_log$estimate__ <- 1 - c_eff$TravDist_log$estimate__
c_eff$TravDist_log$lower__ <- 1 - c_eff$TravDist_log$lower__ 
c_eff$TravDist_log$upper__ <- 1 - c_eff$TravDist_log$upper__

c_eff2$TravDist_log$lower__ <- 1 - c_eff2$TravDist_log$lower__ 
c_eff2$TravDist_log$upper__ <- 1 - c_eff2$TravDist_log$upper__

c_eff3$TravDist_log$lower__ <- 1 - c_eff3$TravDist_log$lower__ 
c_eff3$TravDist_log$upper__ <- 1 - c_eff3$TravDist_log$upper__


bin3 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[3]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Travel time to major cities (log)") +
  geom_ribbon(data = c_eff2$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(bin3)

#### Stunting ####
c_eff$Stunting_log$estimate__ <- 1 - c_eff$Stunting_log$estimate__
c_eff$Stunting_log$lower__ <- 1 - c_eff$Stunting_log$lower__ 
c_eff$Stunting_log$upper__ <- 1 - c_eff$Stunting_log$upper__

c_eff2$Stunting_log$lower__ <- 1 - c_eff2$Stunting_log$lower__ 
c_eff2$Stunting_log$upper__ <- 1 - c_eff2$Stunting_log$upper__

c_eff3$Stunting_log$lower__ <- 1 - c_eff3$Stunting_log$lower__ 
c_eff3$Stunting_log$upper__ <- 1 - c_eff3$Stunting_log$upper__

bin5 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[5]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Prevalence of stunting (log)") +
  geom_ribbon(data = c_eff2$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(bin5)


### Human population density ###
c_eff$PopDens_log$estimate__ <- 1 - c_eff$PopDens_log$estimate__
c_eff$PopDens_log$lower__ <- 1 - c_eff$PopDens_log$lower__ 
c_eff$PopDens_log$upper__ <- 1 - c_eff$PopDens_log$upper__

c_eff2$PopDens_log$lower__ <- 1 - c_eff2$PopDens_log$lower__ 
c_eff2$PopDens_log$upper__ <- 1 - c_eff2$PopDens_log$upper__

c_eff3$PopDens_log$lower__ <- 1 - c_eff3$PopDens_log$lower__ 
c_eff3$PopDens_log$upper__ <- 1 - c_eff3$PopDens_log$upper__


bin6 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[6]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Human population density (log)") +
  geom_ribbon(data = c_eff2$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(bin6)


### Food biomass ###
c_eff$FoodBiomass_log$estimate__ <- 1 - c_eff$FoodBiomass_log$estimate__
c_eff$FoodBiomass_log$lower__ <- 1 - c_eff$FoodBiomass_log$lower__ 
c_eff$FoodBiomass_log$upper__ <- 1 - c_eff$FoodBiomass_log$upper__

c_eff2$FoodBiomass_log$lower__ <- 1 - c_eff2$FoodBiomass_log$lower__ 
c_eff2$FoodBiomass_log$upper__ <- 1 - c_eff2$FoodBiomass_log$upper__

c_eff3$FoodBiomass_log$lower__ <- 1 - c_eff3$FoodBiomass_log$lower__ 
c_eff3$FoodBiomass_log$upper__ <- 1 - c_eff3$FoodBiomass_log$upper__


bin8 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[8]] +
  theme_bw()+ ylab("Probability of being locally extirpated") + xlab("Food biomass (log)") +
  geom_ribbon(data = c_eff2$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(bin8)

### Traded ###
c_eff$Traded$estimate__ <- 1 - c_eff$Traded$estimate__
c_eff$Traded$lower__ <- 1 - c_eff$Traded$lower__ 
c_eff$Traded$upper__ <- 1 - c_eff$Traded$upper__
bin4 <- plot(c_eff, plot=FALSE)[[4]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(bin4)

#### Reserve ####
c_eff$Reserve$estimate__ <- 1 - c_eff$Reserve$estimate__
c_eff$Reserve$lower__ <- 1 - c_eff$Reserve$lower__ 
c_eff$Reserve$upper__ <- 1 - c_eff$Reserve$upper__
bin7 <- plot(c_eff, plot=FALSE)[[7]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(bin7)

### Trophic level ###
c_eff$Trophic_Level$estimate__ <- 1 - c_eff$Trophic_Level$estimate__
c_eff$Trophic_Level$lower__ <- 1 - c_eff$Trophic_Level$lower__ 
c_eff$Trophic_Level$upper__ <- 1 - c_eff$Trophic_Level$upper__
bin9 <- plot(c_eff, plot=FALSE)[[9]] + theme_bw() + ylab("Probability of being locally extirpated")
plot(bin9)

ggpubr::ggarrange(bin1, bin2, bin10, bin3, bin5, bin6, bin8, bin4, bin7, bin9) 
ggsave(plot2, filename="Figures/Model/Bin_Marginal_effects.pdf", width=400, height=313, units="mm", device="pdf")
rm(list=c("mu1", "mu2", "mu3","mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "plot1"))


##### Plotting effects: gaussian #####

c_eff <- conditional_effects(gaussian, method="fitted", prob=0.95)#
c_eff2 <- conditional_effects(gaussian, method="fitted", prob=0.8)
c_eff3 <- conditional_effects(gaussian, method="fitted", prob=0.5)

gau1 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[1]] +
  theme_bw()+ ylab("RR") + xlab("Distance to hunter's access point (log)") +
  geom_ribbon(data = c_eff2$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$DistHunt_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau1)

gau2 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[2]] +
  theme_bw()+ ylab("RR") + xlab("Body mass (log)") +
  geom_ribbon(data = c_eff2$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$BodyMass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau2)

gau9 <- plot(c_eff, plot=FALSE, list(alpha=0.2))[[9]] +theme_bw() + ylab("RR") + 
  xlab("Distance to hunter's access point (log)") + 
  scale_fill_manual(values=c("#A10000", "#D0A277", "#FDAF17"),"Body mass")+
  scale_colour_manual(values=c("#820303", "#BC926B", "#DA9611"), "Body mass") + 
  theme(legend.position=c(0.5,0.9), legend.direction = "horizontal")
plot(gau9)

gau3 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[3]] +
  theme_bw()+ ylab("RR") + xlab("Travel time to major cities (log)") +
  geom_ribbon(data = c_eff2$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$TravDist_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau3)

gau4 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[4]] +
  theme_bw()+ ylab("RR") + xlab("Prevalence of stunting (log)") +
  geom_ribbon(data = c_eff2$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$Stunting_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau4)

gau10 <- plot(c_eff, plot=FALSE, list(alpha=0.2))[[10]] +theme_bw() + ylab("RR") + 
  xlab("Travel time to major cities (log)") + 
  scale_fill_manual(values=c("#A10000", "#D0A277", "#FDAF17"),"Stunting")+
  scale_colour_manual(values=c("#820303", "#BC926B", "#DA9611"), "Stunting") + 
  theme(legend.position=c(0.5,0.05), legend.direction = "horizontal")
plot(gau10)

summary(gaussian)

gau6 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[6]] +
  theme_bw()+ ylab("RR") + xlab("Human population density (log)") +
  geom_ribbon(data = c_eff2$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$PopDens_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau6)

gau7 <- plot(c_eff, plot=FALSE,line_args = list(size=2,alpha=0.1, color="black", fill="#de2d26"))[[7]] +
  theme_bw()+ ylab("RR") + xlab("Food biomass (log)") +
  geom_ribbon(data = c_eff2$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26"))+
  geom_ribbon(data = c_eff3$FoodBiomass_log, aes(ymin=lower__, ymax=upper__,alpha=0.9,fill="#de2d26")) +
  theme(legend.position = "none")
plot(gau7)

gau8 <- plot(c_eff, plot=FALSE)[[8]] + theme_bw() + ylab("RR")
gau5 <- plot(c_eff, plot=FALSE)[[5]] + theme_bw() + ylab("RR")

plot3 <- ggpubr::ggarrange(gau1,gau2, gau9,gau6, gau3, gau4, gau10, gau7, gau8, gau5) 
plot(plot3)
ggsave(plot1, filename="Figures/Model/RR_Marginal_effects.pdf", width=400, height=313, units="mm", device="pdf")
rm(list=c("mu1", "mu2", "mu3","mu4", "mu5", "mu6", "mu7", "mu8", "mu9", "plot1"))



