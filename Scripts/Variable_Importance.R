#Script: estimating varibale importance from hurdle model
# Date: 11th May, 2023

library(brms)
library(ggplot2)

load("Results/hurdle_gam.Rdata")

# estimate R2 (be sure that the dataset is correct. If not, run RR_BRMS.R until line 158)
hurdle_R2<-cor(RR_data$RR, predict(hurdle_gam, newdata = RR_data, summary=TRUE)[,1])^2 #extract "Estimate" column and ^2
hurdle_R2

#create dataframe to save values
out<-data.frame(
  Variables=c('DistHunt_log', 'BodyMass_log', 'TravDist_log', 'Traded', 'Stunting_log', 'PopDens_log', 'Reserve', 'BirdTree_Species', 'Dataset', 'CountryNum'), 
  Name=c('Distance to hunters access point', 'Body mass', 'Travel time to cities', 'Traded', 'Prevalence of stunting', 'Human population density', 'Reserve', 'Species', 'Study', 'Country'), 
  Cat=c("Proxy of hunting pressure","Trait", "Proxy of hunting pressure", "Trait", rep('Proxy of hunting pressure', 3), rep('Random effects', 3)),
  Imp=NA)

for (i in 1:length(out$Variables)){
  print(i)
  var<-as.character(out$Variables[i])
  print(var)
  DX<-RR_data
  if(is.numeric(DX[,var])) { #if numeric
    DX[,var]<-0
    R2<-cor(RR_data$RR, predict(hurdle_gam, newdata = DX, summary=TRUE)[,1])^2
  } else { #if categorical
    if(out$Cat[i]!='Random effects') {
      outCat<-rep(NA, length(unique(RR_data[,var])))
      for (j in 1:length(unique(RR_data[,var]))) {
        DX[,var]<-unique(RR_data[,var])[j]
        outCat[j]<-cor(RR_data$RR, predict(hurdle_gam, newdata = DX, summary=TRUE)[,1])^2
        if(j>10) {break}
      } #close loop
      R2<-mean(outCat, na.rm=TRUE)
    } else { #if a random effect
      DX[,var]<-'NewLevel'
      R2<-cor(RR_data$RR, predict(hurdle_gam, newdata = DX, summary=TRUE, allow_new_levels=TRUE)[,1])^2
    } #close random effect
  } #close categorical
  out[out$Variables==var, 'Imp']<-hurdle_R2-R2
}

out<-out[rev(order(out$Imp)),]
out$ImpPer<-100*out$Imp/sum(out$Imp)
out$ImpPer<-ifelse(out$ImpPer < 0, 0, out$ImpPer)

imp1 <- ggplot(out, aes(y=ImpPer, x=reorder(Name, -ImpPer), fill=Cat)) + geom_col()+ xlab("Predictor") +
  ylab("Predictor contribution on the total variance explained (%)") + labs(fill="Predictor category") +
  scale_fill_manual(values=c("#627ED0","#98BF99", "#F8EEAE")) + theme_bw() +
  theme(legend.position = c(0.9,0.9), legend.background = element_blank()) + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

imp1
#ggsave("Figures/Model/Hurdle_VarImp_RE.pdf", plot=imp1, width=29.7, height=21 , units="cm")

out2<-out[!out$Cat %in% c('Random effects'), ]
out2$ImpPer<-100*out2$Imp/sum(out2$Imp)
out2$ImpPer<-ifelse(out2$ImpPer < 0, 0, out2$ImpPer)

imp2 <- ggplot(out2, aes(y=ImpPer, x=reorder(Name, -ImpPer), fill=Cat)) + geom_col()+ xlab("Predictor") +
  ylab("Predictor contribution on the total variance explained (%)")+ theme_bw()+
  scale_fill_manual(values=c("#627ED0","#F8EEAE")) + labs(fill="Predictor category") +
  theme(legend.position = c(0.9,0.95), legend.background = element_blank()) + ylim(0,55)+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

imp2
#ggsave("Figures/Model/Hurdle_VarImp.pdf", plot=imp2, width=29.7, height=21 , units="cm")

rm(list=c("out","out2", "hurdle_R2", "i", "j", "R2", "var","outCat","DX", "imp2","imp1"))

