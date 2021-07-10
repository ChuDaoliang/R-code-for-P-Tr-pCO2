
##packages#########################################

library(propagate)# for error propagate of pCO2
library(fANCOVA)      # for calculating best smooth for loess
library(RColorBrewer) # for color 
library(ggplot2)      # for data plotting
library(cowplot)      # for multi-graph option

set.seed(22) 

##########################data##############################################################################

#original Î´13C of marine carbonate data
marine.ci.carb<-read.csv("E:statistics/Rdata/pCO2 project/marine ci carb.csv", header = T) 
#original Î´13C of plant data
terrestrial.plant.ci<-read.csv("E:statistics/Rdata/pCO2 project/terrestrial.plant.ci.csv", header = T) 
#original surface sea temperature data
sst<-read.csv("E:statistics/Rdata/pCO2 project/sst.csv", header = T) 


########LOESS fit#########################################################################################
# =============================================================================
# LOESS curves with 0.002 Myr spacing were fitted to Î´13Ccarb, Î´13Cdic,  Î´13C in C3 plants and sea surface temperature.
# At each 0.002 Myr time step, the probability maximum value are identified to be served as input parameters in calculations.
# =============================================================================
#best smooth for LOESS fitting of global marine Î´13Ccarb data
model.best<-loess.as( marine.ci.carb$age, marine.ci.carb$ci ,degree = 2, criterion = "aicc") 
span.marine<-as.numeric(model.best$pars[1])  
terrestrial.span<-0.3

#LOEESS fitting of Î´13C of marine carbonate
age.upper<-251.40
age.lower<-252.104 
spacing<- -0.002
marine.loess<-loess(marine.ci.carb$ci~marine.ci.carb$age, span=span.marine) #loess model
marine.predict<-predict(marine.loess, se=T, newdata=seq(age.lower, age.upper, spacing)) #æ ¹æ®åŽŸæœ‰çš„ageï¼Œåœ¨loessmodel é‡Œpredict sst


#input parameter of terrestrial Î´13Corg 
maxweight<-15
spacing<- -0.002
terrestrial.plant.ci<-terrestrial.plant.ci[order(terrestrial.plant.ci$age, decreasing = T),]

terrestrial.plant.ci<-data.frame(terrestrial.plant.ci, weight= c(rep(1, 107 ), maxweight, rep(1, 49), maxweight, rep(1, 12)))
terrestrial.predict<-predict(loess(terrestrial.plant.ci$ci~terrestrial.plant.ci$age, span=terrestrial.span, weights =terrestrial.plant.ci$weight ), se=T, newdata = seq(age.lower, age.upper, spacing))
terrestrial.plant.ci.t.sd<-terrestrial.predict$se.fit  

#input parameter of Î´13C in marine DIC
marine.dic.ci<-marine.ci.carb
marine.dic.ci$ci<-marine.ci.carb$ci-1 #the fractionation between Î´13C in carbonate and marine DIC is 1â€? (Romanek et al., 1992)  
marine.loess<-loess(marine.dic.ci$ci~marine.dic.ci$age, span=span.marine) 
marine.predict2<-predict(marine.loess, se=T, newdata=seq(age.lower, age.upper, spacing)) 
marine.dic.ci.t.sd<-marine.predict2$se.fit  

#input parameter of Sea surface temperature(sst) 
sst.predict<-predict(loess(sst$sst~sst$sst.age,  span=0.50),se=T, newdata=seq(age.lower, age.upper, spacing))
sst.t.sd<-sst.predict$se.fit  #0  

#data frame of all input parameters
pco2model<-data.frame(age=seq(age.lower, age.upper, spacing), marine.carb.ci=marine.predict$fit, marine.dic.ci=marine.predict2$fit, sst=sst.predict$fit, terrestrial.plant.ci=terrestrial.predict$fit, 
                      marine.carb.ci.sd=marine.predict$se.fit, marine.dic.ci.sd=marine.dic.ci.t.sd, sst.sd=sst.t.sd, terrestrial.plant.ci.sd=terrestrial.plant.ci.t.sd) 



#####################pCO2 calculation and error propagate################################################
# =============================================================================
# t=0, calculation and error propagate of Î´13Ca, âˆ?13C, âˆ?(âˆ?13C), pCO2
# =============================================================================

#fractionation.t0
nsim<-10000
marine.dic.ci.t0.mean<-mean(marine.dic.ci$ci[which(marine.dic.ci$age>252.104)])  #Î´13CDIC(t=0),#the fractionation between Î´13C in carbonate and marine DIC is 1â€? (Romanek et al., 1992)
marine.dic.ci.t0.sd<-sd(marine.dic.ci$ci[which(marine.dic.ci$age>252.104)])


marine.dic.ci.t0<- c(marine.dic.ci.t0.mean, marine.dic.ci.t0.sd)
terrestrial.plant.ci.t0.mean<- -24.42
terrestrial.plant.ci.t0.sd<- 0.5
terrestrial.plant.ci.t0<- c(terrestrial.plant.ci.t0.mean, terrestrial.plant.ci.t0.sd)
sst.t0.mean<-25 #21.9
sst.t0.sd<-0
sst.t0<-c(sst.t0.mean, sst.t0.sd)

nsim<-10000
data1<-cbind(marine.dic.ci.t0, terrestrial.plant.ci.t0, sst.t0)
fractionation.t0.expr<- expression((marine.dic.ci.t0-(0.91*(-0.1141*sst.t0+10.78)+0.08*(-0.052*sst.t0+7.22))-terrestrial.plant.ci.t0)/(1+terrestrial.plant.ci.t0/1000))
res.t0<-propagate(expr=fractionation.t0.expr, data=data1, type="stat", do.sim = T, nsim=nsim)
fractionation.t0<-res.t0$resSIM


# =============================================================================
# t=t, calculation and error propagate of Î´13Ca, âˆ?13C, âˆ?(âˆ?13C), pCO2
# =============================================================================

#Î´13C of atmospheric CO2 (Î´13Ca)
nn<-nrow(pco2model)
atmosphereco2.ci.t.data<-matrix(NA,nrow=nsim,ncol=nn)
i=1
for(i in 1:nn){
  atmosphereco2.ci.t.expr<-expression(marine.dic.ci.t-(0.91*(-0.1141*sst.t+10.78)+0.08*(-0.052*sst.t+7.22)))
  marine.dic.ci.t<-c(pco2model$marine.dic.ci[i], pco2model$marine.dic.ci.sd[i])
  sst.t<-c(pco2model$sst[i], pco2model$sst.sd[i])
  data.atmosphereco2.ci.t<-cbind(marine.dic.ci.t, sst.t)
  res.atmosphereco2.ci.t<-propagate(expr=atmosphereco2.ci.t.expr, data=data.atmosphereco2.ci.t, type="stat", do.sim = T, nsim=nsim)
  atmosphereco2.ci.t.data[,i]<-res.atmosphereco2.ci.t$resSIM
  }

#carbon fractionation between Î´13C of atmospheric CO2 and C3 plant (âˆ?13C)
fractionation.t.data<-matrix(NA,nrow=nsim,ncol=nn)
i=1
for(i in 1:nn){
  fractionation.t.expr<-expression((atmosphereco2.ci.t-terrestrial.plant.ci.t)/(1+terrestrial.plant.ci.t/1000))
  atmosphereco2.ci.t<-atmosphereco2.ci.t.data[,i]
  terrestrial.plant.ci.t<-rnorm(nsim, pco2model$terrestrial.plant.ci[i], pco2model$terrestrial.plant.ci.sd[i])
  data.t<-cbind(atmosphereco2.ci.t, terrestrial.plant.ci.t)
  res.t<-propagate(expr=fractionation.t.expr, data=data.t, type="raw", do.sim = T, nsim=nsim)
  fractionation.t.data[,i]<-res.t$resSIM
  }


#âˆ?13C value between the time of interest (t = t) and a reference time (t = 0) (âˆ?(âˆ?13C))
fractionation.change.expr<- expression(fractionation.t-fractionation.t0)
fractionation.change.data<-matrix(NA,nrow=nsim,ncol=nn)
for(i in 1:nn){
  
  fractionation.t0<-fractionation.t0
  fractionation.t<-fractionation.t.data[,i]
  data.fractionation.change<-cbind(fractionation.t,fractionation.t0 )
  res.fractionation.change<-propagate(expr=fractionation.change.expr, data=data.fractionation.change, type="raw", do.sim = T, nsim=nsim)
  fractionation.change.data[,i]<-res.fractionation.change$resSIM
}

#C parameter
C.expr<-expression((A*4.4)/((A-4.4)*B))
A<-c(28.26,0)  
B<-c(0.223,0.028) 
data.C<-cbind(A,B)
res.C<-propagate(expr=C.expr,data=data.C,type="stat",do.sim = T,nsim=nsim)
C<-res.C$resSIM

#pCO2t
pco2.t.expr<- expression((E*A^2+E*A*B*D+2*E*A*B*C+E*B^2*C*D+E*B^2*C^2+A^2*B*D)/(A^2*B-E*A*B-E*B^2*D-E*B^2*C))
A<-rnorm(nsim,28.26,0)  
B<-rnorm(nsim,0.223,0.028)
C<-C
D<-rnorm(nsim, 425,68) 
pco2.t.data<-matrix(NA,nrow=nsim,ncol=nn)
i=1
for(i in 1:nn){
  E<-fractionation.change.data[,i]
  data.pco2.t<-cbind(A, B, C, D,E)
  res.pco2.t<-propagate(expr=pco2.t.expr, data=data.pco2.t, type="raw", do.sim = T, nsim=nsim)
  pco2.t.data[,i]<-res.pco2.t$resSIM
  
}

# =============================================================================
# function for calculating median, 16th quantile and 84th quantile of 10,000 values for each Î´13Ca, pCO2(t), âˆ?13C and âˆ?(âˆ?13C)
# The invalid pCO2(t) values (i.e., pCO2(t) < 0 or >1000,000 ppm) were excluded.
# =============================================================================

#function for calculating median, 16th quantile and 84th quantile
QUANTILE1<-function(x){
  pco2.t.median<-median(x)
  pco2.t.16th<-quantile(x, 0.16)
  pco2.t.84th<-quantile(x, 0.84)
  return(list(pco2.t.median, pco2.t.16th, pco2.t.84th))
}

#function for calculating median, 16th quantile and 84th quantile without invalid pCO2(t) values < 0 or >1000,000
QUANTILE2<-function(x){
  pco2.t.median<-median(x[which(x>0&x<1000000)])
  pco2.t.16th<-quantile(x[which(x>0&x<1000000)], 0.16)
  pco2.t.84th<-quantile(x[which(x>0&x<1000000)], 0.84)
  return(list(pco2.t.median, pco2.t.16th, pco2.t.84th))
}

#Î´13C of atmospheric CO2 (Î´13Ca)
result.atmosphereco2.ci.t<-apply(atmosphereco2.ci.t.data, 2, QUANTILE1)
result.atmosphereco2.ci.t<-data.frame(matrix(unlist(result.atmosphereco2.ci.t), ncol=3, byrow=T),stringsAsFactors=FALSE)

#carbon fractionation between Î´13C of atmospheric CO2 and C3 plant (âˆ?13C)
result.fractionation.t<-apply(fractionation.t.data, 2, QUANTILE1)
result.fractionation.t<-data.frame(matrix(unlist(result.fractionation.t), ncol=3, byrow=T),stringsAsFactors=FALSE)

#âˆ?13C value between the time of interest (t = t) and a reference time (t = 0) (âˆ?(âˆ?13C))
result.fractionation.change<-apply(fractionation.change.data, 2, QUANTILE1)
result.fractionation.change<-data.frame(matrix(unlist(result.fractionation.change), ncol=3, byrow=T),stringsAsFactors=FALSE)

#pCO2(t)
result.pco2.t<-apply(pco2.t.data, 2, QUANTILE2)
result.pco2.t <- data.frame(matrix(unlist(result.pco2.t), ncol=3, byrow=T),stringsAsFactors=FALSE)

#output data
output.pco2model<-data.frame(  pco2model,
                               atmosphereco2.ci.median=result.atmosphereco2.ci.t$X1, atmosphereco2.ci.16th=result.atmosphereco2.ci.t$X2, atmosphereco2.ci.84th=result.atmosphereco2.ci.t$X3,
                               F.median=result.fractionation.t$X1, F.16th=result.fractionation.t$X2, F.84th=result.fractionation.t$X3, 
                               F.change.median=result.fractionation.change$X1, F.change.16th=result.fractionation.change$X2, F.change.84th=result.fractionation.change$X3, 
                               pco2.median=result.pco2.t$X1, pco2.16th=result.pco2.t$X2, pco2.84th=result.pco2.t$X3
)
write.table(output.pco2model, "C:/Users/Revo/Desktop/output.pco2mode.csv" , col.names=T, row.names=F, sep="," )

output.pco2model<-read.csv("C:/Users/Revo/Desktop/output.pco2mode.csv" , header=T)

