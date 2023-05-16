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
diag <- ggpubr::ggarrange(p1,p2,p3,p4)

ggsave("Figures/Model/Hurdle_Diagnosis.pdf", diag, width = 20, height = 20) 
rm(list=c("p1","p2","p3","p4", "diag") ) #clean environment

# residuals vs fitted values
RR_data2 <- hurdle_gam$data
RR_data2$Residuals <- residuals(hurdle_gam)[,1]
RR_data2$Fitted <- fitted(hurdle_gam)[,1]
plot(RR_data2$Fitted, RR_data2$Residuals)
rm(RR_data2)

#marginal effects: changes in abundance

p1 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("BodyMass_log"), line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(angle=90)) + ylab("RR") + ylim(-2,2) 
p2 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("BodyMass_log"), line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(angle=90)) + ylab("RR") + ylim(-2,2) 
p3 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("BodyMass_log"), line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_text(angle=90)) + ylab("RR") + ylim(-2,2) 

mu1 <- p1 + geom_ribbon(data=p2$data,aes(x=p2$data$x, y=p2$data$predicted, ymin=p2$data$conf.low, ymax=p2$data$conf.high, alpha=0.2)) +
           geom_ribbon(data=p3$data,aes(x=p3$data$x, y=p3$data$predicted, ymin=p3$data$conf.low, ymax=p3$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p4 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("DistHunt_log"), line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) + ylim(-2,2) 
p5 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("DistHunt_log"), line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p6 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("DistHunt_log"), line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) + ylim(-2,2) 

mu2 <- p4 + geom_ribbon(data=p5$data,aes(x=p5$data$x, y=p5$data$predicted, ymin=p5$data$conf.low, ymax=p5$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p6$data,aes(x=p6$data$x, y=p6$data$predicted, ymin=p6$data$conf.low, ymax=p6$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p7 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("TravDist_log"), line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p8 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("TravDist_log"), line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p9 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("TravDist_log"), line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 

mu3 <- p7 + geom_ribbon(data=p8$data,aes(x=p8$data$x, y=p8$data$predicted, ymin=p8$data$conf.low, ymax=p8$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p9$data,aes(x=p9$data$x, y=p9$data$predicted, ymin=p9$data$conf.low, ymax=p9$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p10 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Stunting_log"), line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) + ylim(-2,2) 
p11 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Stunting_log"), line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p12 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Stunting_log"), line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) + ylim(-2,2) 

mu4 <- p10 + geom_ribbon(data=p11$data,aes(x=p11$data$x, y=p11$data$predicted, ymin=p11$data$conf.low, ymax=p11$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p12$data,aes(x=p12$data$x, y=p12$data$predicted, ymin=p12$data$conf.low, ymax=p12$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p13 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("PopDens_log"), line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p14 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("PopDens_log"), line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank())  + ylim(-2,2) 
p15 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("PopDens_log"), line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.title.y = element_blank()) + ylim(-2,2) 

mu5 <- p13 + geom_ribbon(data=p14$data,aes(x=p14$data$x, y=p14$data$predicted, ymin=p14$data$conf.low, ymax=p14$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p15$data,aes(x=p15$data$x, y=p15$data$predicted, ymin=p15$data$conf.low, ymax=p15$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p16 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5, line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)
p17 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5, line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)
p18 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Reserve"), dot.size=5, line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)

mu6 <- p16 + geom_linerange(data=p17$data,aes(x=p17$data$x, y=p17$data$predicted, ymin=p17$data$conf.low, ymax=p17$data$conf.high, alpha=0.2, size=2, color="red")) +
             geom_linerange(data=p18$data,aes(x=p18$data$x, y=p18$data$predicted, ymin=p18$data$conf.low, ymax=p18$data$conf.high, alpha=0.2, size=2), color="firebrick")+ theme(legend.position = "none")


p19 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Traded"), dot.size=5, line.size = 2, ci.lvl=0.5) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)
p20 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Traded"), dot.size=5, line.size = 2, ci.lvl=0.8) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)
p21 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="mu", terms=c("Traded"), dot.size=5, line.size = 2, ci.lvl=0.95) + theme_bw() + theme(plot.title = element_blank(), axis.title.x = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank(),axis.text.x=element_blank())+ ylim(-2,2)

mu7 <- p19 + geom_linerange(data=p20$data,aes(x=p20$data$x, y=p20$data$predicted, ymin=p20$data$conf.low, ymax=p20$data$conf.high, alpha=0.2, size=2, color="red")) +
  geom_linerange(data=p21$data,aes(x=p21$data$x, y=p21$data$predicted, ymin=p21$data$conf.low, ymax=p21$data$conf.high, alpha=0.2, size=2), color="firebrick") + theme(legend.position = "none")

mu <- ggarrange(mu1,mu2,mu3,mu4,mu5, mu6,mu7, nrow = 1, ncol=7) 
mu


#marginal effects:  probability of local extirpationç

p22 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("BodyMass_log"), line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(), axis.text.y=element_text(angle=90)) + ylab("Probability of local extirpation") + ylim(0,1) + xlab("Body mass")
p23<- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("BodyMass_log"), line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(), axis.text.y=element_text(angle=90)) + ylab("Probability of local extirpation") + ylim(0,1) + xlab("Body mass")
p24 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("BodyMass_log"), line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(), axis.text.y=element_text(angle=90)) + ylab("Probability of local extirpation") + ylim(0,1) + xlab("Body mass")

hu1 <- p22 + geom_ribbon(data=p23$data,aes(x=p23$data$x, y=p23$data$predicted, ymin=p23$data$conf.low, ymax=p23$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p24$data,aes(x=p24$data$x, y=p24$data$predicted, ymin=p24$data$conf.low, ymax=p24$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p25 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("DistHunt_log"), line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Distance to hunter's access")
p26 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("DistHunt_log"), line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Distance to hunter's access")
p27 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("DistHunt_log"), line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Distance to hunter's access")

hu2 <- p25 + geom_ribbon(data=p26$data,aes(x=p26$data$x, y=p26$data$predicted, ymin=p26$data$conf.low, ymax=p26$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p27$data,aes(x=p27$data$x, y=p27$data$predicted, ymin=p27$data$conf.low, ymax=p27$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p28 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("TravDist_log"), line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Travel time to cities")
p29 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("TravDist_log"), line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Travel time to cities")
p30 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("TravDist_log"), line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Travel time to cities")

hu3 <- p28 + geom_ribbon(data=p29$data,aes(x=p29$data$x, y=p29$data$predicted, ymin=p29$data$conf.low, ymax=p29$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p30$data,aes(x=p30$data$x, y=p30$data$predicted, ymin=p30$data$conf.low, ymax=p30$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p31 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Stunting_log"), line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Prevalence of stunting")
p32 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Stunting_log"), line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Prevalence of stunting")
p33 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Stunting_log"), line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Prevalence of stunting")

hu4 <- p31 + geom_ribbon(data=p32$data,aes(x=p32$data$x, y=p32$data$predicted, ymin=p32$data$conf.low, ymax=p32$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p33$data,aes(x=p33$data$x, y=p33$data$predicted, ymin=p33$data$conf.low, ymax=p33$data$conf.high, alpha=0.2)) + theme(legend.position = "none")

p34 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("PopDens_log"), line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Human population density")
p35 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("PopDens_log"), line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Human population density")
p36 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("PopDens_log"), line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Human population density")

hu5 <- p34 + geom_ribbon(data=p35$data,aes(x=p35$data$x, y=p35$data$predicted, ymin=p35$data$conf.low, ymax=p35$data$conf.high, alpha=0.2)) +
  geom_ribbon(data=p36$data,aes(x=p36$data$x, y=p36$data$predicted, ymin=p36$data$conf.low, ymax=p36$data$conf.high, alpha=0.2)) + theme(legend.position = "none")


p37 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Reserve"),dot.size=5, line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Reserve")
p38 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Reserve"), dot.size=5,line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Reserve")
p39 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Reserve"),dot.size=5, line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Reserve")

hu6 <- p37 + geom_linerange(data=p38$data,aes(x=p38$data$x, y=p38$data$predicted, ymin=p38$data$conf.low, ymax=p38$data$conf.high, alpha=0.2, size=2, color="red")) +
  geom_linerange(data=p39$data,aes(x=p39$data$x, y=p39$data$predicted, ymin=p39$data$conf.low, ymax=p39$data$conf.high, alpha=0.2, size=2), color="firebrick")+ theme(legend.position = "none")

p40 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Traded"),dot.size=5, line.size = 2, ci.lvl = 0.5) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Traded")
p41 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Traded"), dot.size=5,line.size = 2, ci.lvl = 0.8) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Traded")
p42 <- sjPlot::plot_model(hurdle_gam, type="pred",dpar="hu", terms=c("Traded"),dot.size=5, line.size = 2, ci.lvl = 0.95) + theme_bw() + theme(plot.title = element_blank(),axis.title.y = element_blank(), axis.text.y=element_blank()) + ylim(0,1) + xlab("Traded")

hu7 <- p40 + geom_linerange(data=p41$data,aes(x=p41$data$x, y=p41$data$predicted, ymin=p41$data$conf.low, ymax=p41$data$conf.high, alpha=0.2, size=2, color="red")) +
  geom_linerange(data=p42$data,aes(x=p42$data$x, y=p42$data$predicted, ymin=p42$data$conf.low, ymax=p42$data$conf.high, alpha=0.2, size=2), color="firebrick")+ theme(legend.position = "none")

hu<-ggarrange(hu1,hu2,hu3,hu4,hu5,hu6,hu7, nrow = 1, ncol=7) 

m_eff <- ggarrange(mu,hu, ncol=1)

ggsave("Figures/Model/Hurdle_MarginalEffects.pdf", plot=m_eff, width=34, height=13, units="cm")
rm(list=c(ls(pattern="^p"),ls(pattern="^mu"), ls(pattern="^hu"), "m_eff")) #remove plots p1-p42

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
