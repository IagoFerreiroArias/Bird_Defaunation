#Script: Model diagnosis and marginal effects
#Author: Iago Ferreiro Arias
#Date: 16th March, 2023

library(brms)
library(sjPlot)
library(ggplot2)
library(ggpubr)

load("Results/hurdle_model.RData")
summary(hurdle_model)

#### Check models performance #####
p1<- brms::pp_check(hurdle_model, ndraws=1000) + xlim(0,15)#1
p2<- brms::pp_check(hurdle_model, type="scatter_avg",ndraws=1000) #predicted vs observed
p3<- brms::pp_check(hurdle_model, type="stat_2d", ndraws=1000) # mean and sd of posterior predictions vs observed
p4<- brms::pp_check(hurdle_model, type="loo_pit", ndraws=1000) 
diag <- ggpubr::ggarrange(p1,p2,p3,p4)
diag
rm(list=c("p1","p2","p3","p4"))
performance::r2_bayes(hurdle_model) # Estimate R2

### Estimate NPP and Body Mass values for categories ###
summary(hurdle_model$data$NPP_log)
df_rescale <- read.csv("Results/Rescaling_values.csv") # import values for rescaling predictors
(log10(4.5) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #3kg = 1.49
(log10(2.5) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #2kg = 1.37
(log10(1.5) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #700gr = 0.71
(log10(0.2) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #200gr = -0.278
(log10(0.015) - df_rescale[which(df_rescale$Predictor=="BodyMass"),"Mean"])/df_rescale[which(df_rescale$Predictor=="BodyMass"),"SD"] #30gr = -1.26
rm(df_rescale) # Create interaction plots with three different values of body mass: 3kg (1.46), 700 gr (0.71) and 30 gr (-0.90)
quantile(hurdle_model$data$NPP_log, probs=seq(0,1,0.25), na.rm=TRUE)


######### Changes in abundances ###############

mu1 <- sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("Stunting_log", "BodyMass_log[1.49, -0.278, -1.35]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) + labs(colour="Body Mass") +
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+
  scale_fill_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ xlab("Prevalence of stunting") + 
  theme(plot.title = element_blank(), legend.position ="None") + geom_hline(yintercept=0, linetype="dashed")

mu2 <- sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("TravDist_log", "BodyMass_log[1.49, -0.278, -1.35]]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) + labs(colour="Body Mass") +
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+
  scale_fill_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09")) + xlab("Travel time to major cities")+
  theme(plot.title = element_blank(), legend.position ="None") + geom_hline(yintercept=0, linetype="dashed")

mu3 <- sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("DistHunt_log", "BodyMass_log[1.49, -0.278, -1.35]]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank(),legend.position ="None") + ylab("log(RR)") + ylim(-2,2) + labs(colour="Body Mass") +
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+
  scale_fill_manual(labels = c("0.03kg","0.7kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09")) + xlab("Distance to hunter access points")+
  theme(plot.title = element_blank()) + geom_hline(yintercept=0, linetype="dashed")

mu4 <-sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("DistHunt_log", "NPP_log[0.77, 0, -0.7]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) +  xlab("Distance to hunter access points") +
  theme(plot.title = element_blank(),legend.position ="None") + geom_hline(yintercept=0, linetype="dashed") + labs(colour="NPP") +
  scale_color_manual(labels = c("Low","Moderate", "High"), values=c("#F4D166", "#95B958", "#38884C"))+
  scale_fill_manual(labels = c("Low","Moderate", "High"), values=c("#F4D166", "#95B958", "#38884C"))

mu5<- sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("PopDens_log"), line.size = 1.5, dot.size = 4) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) +  xlab("Human Population Density")+
  theme(plot.title = element_blank(), legend.position ="None") + geom_hline(yintercept=0, linetype="dashed")

mu6 <- sjPlot::plot_model(hurdle_model, type="pred", dpar="mu", terms=c("Reserve"), line.size = 1.5, dot.size = 4) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("log(RR)") + ylim(-2,2) +  
  theme(plot.title = element_blank(), legend.position ="None") + geom_hline(yintercept=0, linetype="dashed")

mu<-ggpubr::ggarrange( mu1, mu2,mu3,mu6, mu5, mu4, ncol=3, nrow=2, widths = c(0.6,0.6,0.6,1.2,0.6,0.6,0.6,2))
mu #save 7x10 inches 

##### Probability of extinction ####################

hu1 <- sjPlot::plot_model(hurdle_model, type="pred",dpar="hu",terms=c("Stunting_log[all]", "BodyMass_log[1.49, -0.278, -1.35]]"), line.size = 1.5)  + 
  theme_bw() + theme(plot.title = element_blank(), legend.position ="None") + ylab("Probability of local extinction") + ylim(0,1) +labs(colour="Body Mass") +
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ 
  scale_fill_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ xlab("Prevalence of stunting") 

hu2 <- sjPlot::plot_model(hurdle_model, type="pred",dpar="hu",terms=c("TravDist_log[all]", "BodyMass_log[1.49, -0.278, -1.35]]"), line.size = 1.5)  + 
  theme_bw() + theme(plot.title = element_blank(), legend.position ="None") + ylab("Probability of local extinction") + ylim(0,1) +labs(colour="Body Mass") +
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ 
  scale_fill_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ xlab("Travel time to major cities") 

hu3<-sjPlot::plot_model(hurdle_model, type="pred", dpar="hu", terms=c("DistHunt_log[all]", "BodyMass_log[1.49, -0.278, -1.35]]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank(),legend.position ="None") + ylab("Probability of local extinction") + ylim(0,1)+
  scale_color_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09"))+ labs(colour="Body Mass") +
  scale_fill_manual(labels = c("0.03kg","0.2kg", "3kg"), values=c("#3B99B1", "#EAC728", "#F13D09")) + xlab("Distance to hunter access points")

hu4<-sjPlot::plot_model(hurdle_model, type="pred", dpar="hu", terms=c("DistHunt_log[all]", "NPP_log[0.77, 0, -0.7]"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extinction") + ylim(0,1) + 
  theme(plot.title = element_blank(),legend.position ="None") + xlab("Distance to hunters access point") +labs(colour="NPP") +
  scale_color_manual(labels = c("Low","Moderate", "High"), values=c("#F4D166", "#95B958", "#38884C"))+
  scale_fill_manual(labels = c("Low","Moderate", "High"), values=c("#F4D166", "#95B958", "#38884C"))

hu5<-sjPlot::plot_model(hurdle_model, type="pred", dpar="hu", terms=c("PopDens_log"), line.size = 1.5) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extinction") + ylim(0,1) +
  theme(plot.title = element_blank()) + xlab("Human Population Density") 

hu6<-sjPlot::plot_model(hurdle_model, type="pred", dpar="hu", terms=c("Reserve"), line.size = 1.5, dot.size=4) + 
  theme_bw() + theme(plot.title = element_blank()) + ylab("Probability of local extinction") + ylim(0,1) +
  xlab("Reserve") + theme(plot.title = element_blank(),legend.position ="None" )

hu<-ggpubr::ggarrange(hu1, hu2,hu3,hu6, hu5, hu4, ncol=3, nrow=2)
hu #probability of extirpation

#Calculate probability of direction of effect
bayestestR::p_direction(hurdle_model, method="direct")

#Extract summary of the model
summary(hurdle_model_out)

#Clean environment
rm(list=c("mu1", "mu2","mu3","mu4","mu5","mu6","mu7", "mu8", "mu9","mu",
          "hu1", "hu2","hu3","hu4","hu5","hu6","hu7","hu8", "hu9","hu"))