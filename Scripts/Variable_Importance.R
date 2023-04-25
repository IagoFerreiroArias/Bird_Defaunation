##### ESTIMATING VARIABLE IMPORTANCE ####

load("Results/binomial_model.Rdata") #Binomial model (probability of local extirpation)

bin_R2<-cor(RR_data$bin, predict(binomial, newdata = RR_data, summary=TRUE)[,1])^2 #extract "Estimate" column and ^2
bin_R2
#fullR2<-cor(DATA$Dens, predict(mod, newdata = DATA, re_formula=NA, summary=TRUE)[,1])^2

out<-data.frame(
  Variables=c('DistHunt_log', 'BodyMass_log', 'TravDist_log', 'Traded', 'Stunting_log', 'PopDens_log', 'Reserve', 'FoodBiomass_log', 'Trophic_Level', 'BirdTree_Species', 'Dataset', 'CountryNum'), 
  Name=c('Distance hunters access point', 'Body mass', 'Travel time to major cities', 'Traded', 'Prevalence of stunting', 'Human population density', 'Reserve', 'Food biomass', 'Trophic_Level', 'Species', 'Study', 'Country'), 
  Cat=c("Proxy of hunting pressure","Trait", "Proxy of hunting pressure", "Trait", rep('Proxy of hunting pressure', 4), "Trait", rep('Random effects', 3)),
  Imp=NA)

for (i in 1:length(out$Variables)){
  print(i)
  var<-as.character(out$Variables[i])
  print(var)
  DX<-RR_data
  if(is.numeric(DX[,var])) { #if numeric
    DX[,var]<-0
    R2<-cor(RR_data$bin, predict(binomial, newdata = RR_data, summary=TRUE)[,1])^2
  } else { #if categorical
    if(out$Cat[i]!='Random effects') {
      outCat<-rep(NA, length(unique(RR_data[,var])))
      for (j in 1:length(unique(RR_data[,var]))) {
        DX[,var]<-unique(RR_data[,var])[j]
        outCat[j]<-cor(RR_data$bin, predict(binomial, newdata = DX, summary=TRUE)[,1])^2
        if(j>10) {break}
      } #close loop
      R2<-mean(outCat, na.rm=TRUE)
    } else { #if a random effect
      DX[,var]<-'NewLevel'
      R2<-cor(RR_data$bin, predict(binomial, newdata = DX, summary=TRUE, allow_new_levels=TRUE)[,1])^2
    } #close random effect
  } #close categorical
  out[out$Variables==var, 'Imp']<-bin_R2-R2
}

out<-out[rev(order(out$Imp)),]
out$ImpPer<-100*out$Imp/sum(out$Imp)
out$ImpPer<-ifelse(out$ImpPer < 0, 0, out$ImpPer)

out2<-out[!out$Cat %in% c('Random effects'), ]
out2$ImpPer<-100*out2$Imp/sum(out2$Imp)
out2$ImpPer<-ifelse(out2$ImpPer < 0, 0, out2$ImpPer)

write.csv(out, "Results/Binomial_VarImp1.csv")
write.csv(out2, "Results/Binomial_VarImp2.csv")

load("Results/gaussian_model.Rdata") #gaussian model (changes in abundance)

new_data <- RR_data[which(RR_data$RR_log !='-Inf'),]

gau_R2<-cor(new_data$RR_log, predict(gaussian, newdata = new_data, summary=TRUE)[,1])^2 #extract "Estimate" column and ^2
gau_R2

#fullR2<-cor(DATA$Dens, predict(mod, newdata = DATA, re_formula=NA, summary=TRUE)[,1])^2

formula(gaussian)
out<-data.frame(
  Variables=c('DistHunt_log', 'BodyMass_log', 'TravDist_log', 'Stunting_log', 'Reserve', 'PopDens_log', 'FoodBiomass_log', 'Traded', 'BirdTree_Species', 'Dataset', 'CountryNum'), 
  Name=c('Distance hunters access point', 'Body mass', 'Travel time to major cities', 'Prevalence of stunting', 'Reserve', 'Human population density', 'Food biomass', 'Traded', 'Species', 'Study', 'Country'), 
  Cat=c("Proxy of hunting pressure","Trait", rep('Proxy of hunting pressure', 5), "Trait", rep('Random effects', 3)),
  Imp=NA)

for (i in 1:length(out$Variables)){
  print(i)
  var<-as.character(out$Variables[i])
  print(var)
  DX<-new_data
  if(is.numeric(DX[,var])) { #if numeric
    DX[,var]<-0
    R2<-cor(new_data$RR_log, predict(gaussian, newdata = new_data, summary=TRUE)[,1])^2
  } else { #if categorical
    if(out$Cat[i]!='Random effects') {
      outCat<-rep(NA, length(unique(new_data[,var])))
      for (j in 1:length(unique(new_data[,var]))) {
        DX[,var]<-unique(new_data[,var])[j]
        outCat[j]<-cor(new_data$RR_log, predict(gaussian, newdata = DX, summary=TRUE)[,1])^2
        if(j>10) {break}
      } #close loop
      R2<-mean(outCat, na.rm=TRUE)
    } else { #if a random effect
      DX[,var]<-'NewLevel'
      R2<-cor(new_data$RR_log, predict(gaussian, newdata = DX, summary=TRUE, allow_new_levels=TRUE)[,1])^2
    } #close random effect
  } #close categorical
  out[out$Variables==var, 'Imp']<-gau_R2-R2
}

out<-out[rev(order(out$Imp)),]
out$ImpPer<-100*out$Imp/sum(out$Imp)
out$ImpPer<-ifelse(out$ImpPer < 0, 0, out$ImpPer)

out2<-out[!out$Cat %in% c('Random effects'), ]
out2$ImpPer<-100*out2$Imp/sum(out2$Imp)
out2$ImpPer<-ifelse(out2$ImpPer < 0, 0, out2$ImpPer)


write.csv(out, "Results/Gaussian_VarImp1.csv")
write.csv(out2, "Results/Gaussian_VarImp2.csv")
