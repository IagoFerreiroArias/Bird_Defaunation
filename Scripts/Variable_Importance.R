#Script: estimating varibale importance from hurdle model
# Date: 11th May, 2023

library(brms)
library(ggplot2)
library(dplyr)

load("Results/hurdle_model.Rdata")
RR_data <- hurdle_model$data

# estimate R2 (be sure that the dataset is correct. If not, run RR_BRMS.R until line 158)
hurdle_R2<-cor(RR_data$RR, predict(hurdle_model, newdata = RR_data, summary=TRUE)[,1])^2 #extract "Estimate" column and ^2
hurdle_R2

#create dataframe to save values
out<-data.frame(
  Variables=c('DistHunt_log','BodyMass_log', 'TravDist_log', 'Stunting_log', 'PopDens_log', 'NPP_log','Reserve', 'BirdTree_Species', 'CountryNum'), 
  Name=c('Distance to hunters access point', 'Body mass', 'Travel time to cities','Prevalence of stunting', 'Human population density', 'Net Primary Productivity','Reserve', 'Species', 'Country'), 
  Cat=c('Proxy of hunting pressure',"Trait", rep('Proxy of hunting pressure', 5), rep('Random effects', 2)),
  Imp=NA)

for (i in 1:length(out$Variables)){
  print(i)
  var<-as.character(out$Variables[i])
  print(var)
  DX<-RR_data
  if(is.numeric(DX[,var])) { #if numeric
    DX[,var]<-0
    R2<-cor(RR_data$RR, predict(hurdle_model, newdata = DX, summary=TRUE)[,1])^2
  } else { #if categorical
    if(out$Cat[i]!='Random effects') {
      outCat<-rep(NA, length(unique(RR_data[,var])))
      for (j in 1:length(unique(RR_data[,var]))) {
        DX[,var]<-unique(RR_data[,var])[j]
        outCat[j]<-cor(RR_data$RR, predict(hurdle_model, newdata = DX, summary=TRUE)[,1])^2
        if(j>10) {break}
      } #close loop
      R2<-mean(outCat, na.rm=TRUE)
    } else { #if a random effect
      DX[,var]<-'NewLevel'
      R2<-cor(RR_data$RR, predict(hurdle_model, newdata = DX, summary=TRUE, allow_new_levels=TRUE)[,1])^2
    } #close random effect
  } #close categorical
  out[out$Variables==var, 'Imp']<-hurdle_R2-R2
}

out<-out[rev(order(out$Imp)),]
out$ImpPer<-100*out$Imp/sum(out$Imp)
out$ImpPer<-ifelse(out$ImpPer < 0, 0, out$ImpPer)


imp1 <- ggplot(out, aes(y=ImpPer, x=reorder(Name, ImpPer), fill=Cat)) + geom_col()+ xlab("Predictor") +
  ylab("Predictor contribution on the total variance explained (%)") + labs(fill="Predictor category") +
  scale_fill_manual(values=c("#FF5733","#3989B9", "#FFC300")) + theme_bw() +
  theme(legend.position = c(0.75,0.15), legend.background = element_blank(),
        legend.key.size = unit(1, 'cm')) + coord_flip()

imp1

out <- out %>% dplyr::arrange(desc(ImpPer))
out$ImpPer <- round(out$ImpPer,2)

ggplot(out, aes(x="", y=reorder(Name, ImpPer), fill=Name)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_void() + 
  scale_fill_manual("#3B99B1", "#36A5AA", "#59B29D", "#96BE95",
                    "#D9CB80", "#EAC728", "#EB6100", "#F5191C")
  

#ggsave("Figures/Model/Hurdle_VarImp_RE.pdf", plot=imp1, width=29.7, height=21 , units="cm")

out2<-out[!out$Cat %in% c('Random effects'), ]
out2$ImpPer<-100*out2$Imp/sum(out2$Imp)
out2$ImpPer<-ifelse(out2$ImpPer < 0, 0, out2$ImpPer)

imp2 <- ggplot(out2, aes(y=ImpPer, x=reorder(Name, -ImpPer), fill=Cat)) + geom_col()+ xlab("Predictor") +
  ylab("Predictor contribution on the total variance explained (%)")+ theme_bw()+
  scale_fill_manual(values=c("#FF5733","#FFC300")) + labs(fill="Predictor category") +
  theme(legend.position = c(0.9,0.95), legend.background = element_blank())
  

imp2
#ggsave("Figures/Model/Hurdle_VarImp.pdf", plot=imp2, width=29.7, height=21 , units="cm")

rm(list=c("out","out2", "hurdle_R2", "i", "j", "R2", "var","outCat","DX", "imp2","imp1"))

