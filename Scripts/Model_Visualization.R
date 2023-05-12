<<<<<<< Updated upstream
#Script: Model diagnosis and marginal effects
#Author: Iago Ferreiro Arias
#Date: 16th March, 2023

library(brms)
library(sjPlot)
library(ggplot2)
library(ggpubr)

# import hurdle model
load("Results/hurdle_gam.RData")

# diagnosis plots

p1<- brms::pp_check(hurdle_gam, ndraws=1000) + xlim(0,15)#1
p2<- brms::pp_check(hurdle_gam, type="scatter_avg",ndraws=1000) #predicted vs observed
p3<- brms::pp_check(hurdle_gam, type="stat_2d", ndraws=1000) # mean and sd of posterior predictions vs observed
p4<- brms::pp_check(hurdle_gam, type="loo_pit", ndraws=1000) 
ggpubr::ggarrange(p1,p2,p3,p4)

rm(list=c("p1","p2","p3","p4"))


# residuals vs fitted values
RR_data2 <- hurdle_gam$data
RR_data2$Residuals <- residuals(hurdle_gam)[,1]
RR_data2$Fitted <- fitted(hurdle_gam)[,1]
plot(RR_data2$Fitted, RR_data2$Residuals)
rm(RR_data2)

#marginal effects: changes in abundance

mu1 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("BodyMass_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu2 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("DistHunt_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu3 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu4 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu5 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu6 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu7 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
ggarrange(mu1,mu2,mu3,mu4,mu5, mu6,mu7, nrow = 2, ncol=4) 

rm(list=c("mu1","mu2","mu3", "mu4","mu5","mu6","mu7"))

#marginal effects:  probability of local extirpation
hu1 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("BodyMass_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu2 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("DistHunt_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu3 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu4 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu5 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu6 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("Probability of local extirpation") + ylim(0,1)
hu7 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("Probability of local extirpation") + ylim(0,1)
ggarrange(hu1,hu2,hu3,hu4,hu5,hu6,hu7, nrow = 2, ncol=4) 
rm(list=c("hu1","hu2","hu3", "hu4","hu5","hu6","hu7", "hurdle_gam"))


##### Hurdle gam with interaction between body mass and distance to hunter's access points
load("Results/hurdle_gam_int.RData")

# diagnosis plots
p1<- brms::pp_check(hurdle_gam_int, ndraws=1000) + xlim(0,15)#1
p2<- brms::pp_check(hurdle_gam_int, type="scatter_avg",ndraws=1000) #predicted vs observed
p3<- brms::pp_check(hurdle_gam_int, type="stat_2d", ndraws=1000) # mean and sd of posterior predictions vs observed
p4<- brms::pp_check(hurdle_gam_int, type="loo_pit", ndraws=1000) 
ggpubr::ggarrange(p1,p2,p3,p4)

#warning messages: In smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : basis dimension, k, increased to minimum possible
# interaction between body mass and distance to hunter's access point not supported?

rm(list=c("p1","p2","p3","p4"))


# residuals vs fitted values
RR_data2 <- hurdle_gam_int$data
RR_data2$Residuals <- residuals(hurdle_gam_int)[,1]
RR_data2$Fitted <- fitted(hurdle_gam_int)[,1]
plot(RR_data2$Fitted, RR_data2$Residuals)
rm(RR_data2)


#stablish reference values of body mass in order to plot interactions

df_rescale <- read.csv("Results/Rescaling_values.csv") # import values for rescaling predictors
(log10(3) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #4kg = 1.46
(log10(2.5) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #2kg = 1.37
(log10(0.7) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #700gr = 0.71
(log10(0.03) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #30gr = -0.9
(log10(0.015) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #30gr = -1.26
rm(df_rescale) # Create interaction plots with three different values of body mass: 3kg (1.46), 700 gr (0.71) and 30 gr (-0.90)

# marginal effects: changes in abundance
mu1 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("DistHunt_log", "BodyMass_log[1.46,0.71,-0.9")) + 
               theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) + labs(colour="Body Mass") +
               scale_color_discrete(labels = c("0.03kg","0.7kg", "2.5kg", "4kg"))
mu2 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu3 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu4 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu5 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu6 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu7 <- ggarrange(mu2,mu3,mu4,mu5,mu6) 
ggarrange(mu7, mu1)

# marginal effects: probability of local extirpation
hu1 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("DistHunt_log", "BodyMass_log[1.46,0.71,-0.9")) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1) + labs(colour="Body Mass") +
  scale_color_discrete(labels = c("0.03kg","0.7kg", "2.5kg", "4kg"))
hu2 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1)
hu3 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)")+ ylim(0,1)
hu4 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1)
hu5 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)")+ ylim(0,1)
hu6 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(0,1)
hu7 <- ggarrange(hu2,hu3,hu4,hu5,hu6) 
ggarrange(hu7, hu1)
rm(list=c("hu1","hu2","hu3", "hu4","hu5","hu6","hu7", "hurdle_gam_int"))
# interaction between distance to hunter's access point and body mass is not supported
=======
#Script: Model diagnosis and marginal effects
#Author: Iago Ferreiro Arias
#Date: 16th March, 2023

library(brms)
library(sjPlot)
library(ggplot2)
library(ggpubr)

# import hurdle model
load("Results/hurdle_gam.RData")

# diagnosis plots

p1<- brms::pp_check(hurdle_gam, ndraws=1000) + xlim(0,15)#1
p2<- brms::pp_check(hurdle_gam, type="scatter_avg",ndraws=1000) #predicted vs observed
p3<- brms::pp_check(hurdle_gam, type="stat_2d", ndraws=1000) # mean and sd of posterior predictions vs observed
p4<- brms::pp_check(hurdle_gam, type="loo_pit", ndraws=1000) 
ggpubr::ggarrange(p1,p2,p3,p4)

rm(list=c("p1","p2","p3","p4"))


# residuals vs fitted values
RR_data2 <- hurdle_gam$data
RR_data2$Residuals <- residuals(hurdle_gam)[,1]
RR_data2$Fitted <- fitted(hurdle_gam)[,1]
plot(RR_data2$Fitted, RR_data2$Residuals)
rm(RR_data2)

#marginal effects: changes in abundance

mu1 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("BodyMass_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu2 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("DistHunt_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu3 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu4 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu5 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu6 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu7 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
ggarrange(mu1,mu2,mu3,mu4,mu5, mu6,mu7, nrow = 2, ncol=4) 

rm(list=c("mu1","mu2","mu3", "mu4","mu5","mu6","mu7"))

#marginal effects:  probability of local extirpation
hu1 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("BodyMass_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu2 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("DistHunt_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu3 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu4 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu5 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extirpation") + ylim(0,1)
hu6 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("Probability of local extirpation") + ylim(0,1)
hu7 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("Probability of local extirpation") + ylim(0,1)
ggarrange(hu1,hu2,hu3,hu4,hu5,hu6,hu7, nrow = 2, ncol=4) 
rm(list=c("hu1","hu2","hu3", "hu4","hu5","hu6","hu7", "hurdle_gam"))


##### Hurdle gam with interaction between body mass and distance to hunter's access points
load("Results/hurdle_gam_int.RData")

# diagnosis plots
p1<- brms::pp_check(hurdle_gam_int, ndraws=1000) + xlim(0,15)#1
p2<- brms::pp_check(hurdle_gam_int, type="scatter_avg",ndraws=1000) #predicted vs observed
p3<- brms::pp_check(hurdle_gam_int, type="stat_2d", ndraws=1000) # mean and sd of posterior predictions vs observed
p4<- brms::pp_check(hurdle_gam_int, type="loo_pit", ndraws=1000) 
ggpubr::ggarrange(p1,p2,p3,p4)

#warning messages: In smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : basis dimension, k, increased to minimum possible
# interaction between body mass and distance to hunter's access point not supported?

rm(list=c("p1","p2","p3","p4"))


# residuals vs fitted values
RR_data2 <- hurdle_gam_int$data
RR_data2$Residuals <- residuals(hurdle_gam_int)[,1]
RR_data2$Fitted <- fitted(hurdle_gam_int)[,1]
plot(RR_data2$Fitted, RR_data2$Residuals)
rm(RR_data2)


#stablish reference values of body mass in order to plot interactions

df_rescale <- read.csv("Results/Rescaling_values.csv") # import values for rescaling predictors
(log10(3) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #4kg = 1.46
(log10(2.5) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #2kg = 1.37
(log10(0.7) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #700gr = 0.71
(log10(0.03) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #30gr = -0.9
(log10(0.015) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #30gr = -1.26
rm(df_rescale) # Create interaction plots with three different values of body mass: 3kg (1.46), 700 gr (0.71) and 30 gr (-0.90)

# marginal effects: changes in abundance
mu1 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("DistHunt_log", "BodyMass_log[1.46,0.71,-0.9")) + 
               theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) + labs(colour="Body Mass") +
               scale_color_discrete(labels = c("0.03kg","0.7kg", "2.5kg", "4kg"))
mu2 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu3 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu4 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2)
mu5 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu6 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="mu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(-2,2)
mu7 <- ggarrange(mu2,mu3,mu4,mu5,mu6) 
ggarrange(mu7, mu1)

# marginal effects: probability of local extirpation
hu1 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("DistHunt_log", "BodyMass_log[1.46,0.71,-0.9")) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1) + labs(colour="Body Mass") +
  scale_color_discrete(labels = c("0.03kg","0.7kg", "2.5kg", "4kg"))
hu2 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("TravDist_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1)
hu3 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Stunting_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)")+ ylim(0,1)
hu4 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("PopDens_log")) + theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(0,1)
hu5 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Reserve"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)")+ ylim(0,1)
hu6 <- sjPlot::plot_model(hurdle_gam_int, type="pred",dpar="hu", terms=c("Traded"), dot.size=5) + theme_bw() + theme(plot.title = element_blank())+ ylab("log(RR)") + ylim(0,1)
hu7 <- ggarrange(hu2,hu3,hu4,hu5,hu6) 
ggarrange(hu7, hu1)
rm(list=c("hu1","hu2","hu3", "hu4","hu5","hu6","hu7", "hurdle_gam_int"))
# interaction between distance to hunter's access point and body mass is not supported
>>>>>>> Stashed changes
