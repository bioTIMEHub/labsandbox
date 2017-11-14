# Script for House et al: Moving to 3D: relationships between coral planar area, surface area and volume
# prepared by Jenny House and Maria Dornelas
# 13 September 2016


## DATA MANAGEMENT
data<-(read.table("House_et_al_Data.csv",header=T, row.names=1, sep=","))
data[,1]<-factor(data[,1],ordered=F) #specimen
newdata_means<-aggregate(data[, 2:9], list(data$specimen), mean)
#growthformM is simplified growth form into branching, encrusting and massive
datag<-c()
for(spec in unique(data$specimen))
{
  datag<-c(datag,unique(as.character(data$growthformM[data$specimen==spec])))
}
newdata_means$growthformM=factor(datag)

newdata_means$ctvollog<-log(newdata_means$ct_vol_total)
newdata_means$cttsalog<-log(newdata_means$ct_sa_total)
newdata_means$ctlsalog<-log(newdata_means$ct_sa_live)

newdata_means$planartsalog<-log(newdata_means$planar_total)
newdata_means$planarlsalog<-log(newdata_means$planar_live)

newdata_means$photovollog<-log(newdata_means$photo_vol_total)
newdata_means$phototsalog<-log(newdata_means$photo_sa_total)
newdata_means$photolsalog<-log(newdata_means$photo_sa_live)

str(newdata_means)
View(newdata_means)
attach(newdata_means)

#######################################
# t-tests comparing CT and PH metrics

t.test(cttsalog,phototsalog,paired=TRUE)    ### p=0.002944
t.test(ctvollog,photovollog,paired=TRUE)    ### p=0.0002431
t.test(ctlsalog,photolsalog,paired=TRUE)    ### 4.554e-05


#####################################################################
## MORPHOTYPE ANALYSIS for planar area
# models for total area as a function of planar area

mtsal<-lm(cttsalog~planartsalog)
mtsagfil<-lm(cttsalog~planartsalog+growthformM)
mtsagfsil<-lm(cttsalog~planartsalog*growthformM)

summary(mtsal) 
summary(mtsagfil)
summary(mtsagfsil) 

AIC(mtsal, mtsagfil, mtsagfsil) 
ctsap<-coefficients(mtsal)
########################################################
# models for live area as a function of live planar area

mlsal<-lm(ctlsalog~planarlsalog)
mlsagfil<-lm(ctlsalog~planarlsalog+growthformM)
mlsagfsil<-lm(ctlsalog~planarlsalog*growthformM)

summary(mlsal) 
summary(mlsagfil)
summary(mlsagfsil) 

AIC(mlsal, mlsagfil, mlsagfsil) 
clsap<-coefficients(mlsagfsil)

# best model is full model
##########################################
# models for volume as a function of  planar area

mvoll<-lm(ctvollog~planartsalog)
mvolgfil<-lm(ctvollog~planartsalog+growthformM)
mvolgfsil<-lm(ctvollog~planartsalog*growthformM)

summary(mvoll) 
summary(mvolgfil)
summary(mvolgfsil) 


AIC(mvoll, mvolgfil, mvolgfsil) 
cvsap<-coefficients(mvolgfil)
# best model has intercepts per group (some support for full model)

##########################################

## MORPHOTYPE ANALYSIS for photogammetry
# models for total area as a function of planar area

mtphsal<-lm(cttsalog~phototsalog)
mtphsagfil<-lm(cttsalog~phototsalog+growthformM)
mtphsagfsil<-lm(cttsalog~phototsalog*growthformM)

summary(mtphsal) 
summary(mtphsagfil)
summary(mtphsagfsil) 


AIC(mtphsal, mtphsagfil, mtphsagfsil) 
ctsaph<-coefficients(mtphsal)

########################################################
# models for live area as a function of live photo area

mlphsal<-lm(ctlsalog~photolsalog)
mlphsagfil<-lm(ctlsalog~photolsalog+growthformM)
mlphsagfsil<-lm(ctlsalog~photolsalog*growthformM)

summary(mlphsal) 
summary(mlphsagfil)
summary(mlphsagfsil) 


AIC(mlphsal, mlphsagfil, mlphsagfsil) 
clsaph<-coefficients(mlphsal)
# best model is simplest
##########################################
# models for volume as a function of  photo area

mphvoll<-lm(ctvollog~photovollog)
mphvolgfil<-lm(ctvollog~photovollog+growthformM)
mphvolgfsil<-lm(ctvollog~photovollog*growthformM)

summary(mphvoll) 
summary(mphvolgfil)
summary(mphvolgfsil) 


AIC(mphvoll, mphvolgfil, mphvolgfsil) 
cvsaph<-coefficients(mphvolgfsil)
# best model is full model

##########################################

#############################################################################################################
## Figure 1
colcode<-as.character(growthformM)
colcode[colcode=="branching"]<-"red"
colcode[colcode=="encrusting"]<-"blue4"
colcode[colcode=="massive"]<-"chartreuse3"

par(mfrow=c(3,1))
par(mar=c(4,5,1,2))
par(oma=c(0,0,0,0))


plot(cttsalog~ planartsalog, xlim=c(3, 8),ylim=c(3,9),ylab=expression(paste("Log [CT TSA ", "(mm"^{2},")]")),
     xlab=expression(paste("Log [PL TSA ", "(mm"^{2},")]")),cex.lab=1,yaxs = "i", xaxs="i",bty="n", pch=16, col=colcode)     
abline(0,1) 
abline(ctsap, lty=2)
text(3.5, 8.5, labels="A")


plot(ctlsalog~ planarlsalog, xlim=c(3, 8),ylim=c(3,9), ylab=expression(paste("Log [CT LSA ", "(mm"^{2},")]")),
     xlab=expression(paste("Log [PL LSA ", "(mm"^{2},")]")),cex.lab=1,yaxs = "i", xaxs="i",bty="n", pch=16, col=colcode)     
abline(0,1) 
abline(c(clsap[1],clsap[2]), lty=2, col="red")
abline(c(clsap[1]+clsap[3],clsap[2]+clsap[5]), lty=2, col="blue4")
abline(c(clsap[1]+clsap[4],clsap[2]+clsap[6]), lty=2, col="chartreuse3")
text(3.5, 8.5, labels="B")

plot(ctvollog~ planartsalog, xlim=c(2, 8),ylim=c(2,8), ylab=expression(paste("Log [CT Vol ", "(mm"^{3},")]")),
     xlab=expression(paste("Log [PL TSA ","(mm"^{2},")]")),cex.lab=1,yaxs = "i", xaxs="i",bty="n", pch=16, col=colcode)     
abline(0,1.5)
abline(c(cvsap[1],cvsap[2]), lty=2, col="red")
abline(c(cvsap[1]+cvsap[3],cvsap[2]), lty=2, col="blue4")
abline(c(cvsap[1]+cvsap[4],cvsap[2]), lty=2, col="chartreuse3")
text(2.5, 7, labels="C")



###########################################################################
# Figure 2
plot(cttsalog~ phototsalog, xlim=c(4, 9),ylim=c(4,9), 
     xlab=expression(paste("Log [PH TSA ", "(mm"^{2},")]")),
     ylab=expression(paste("Log [CT TSA ", "(mm"^{2},")]")),cex.lab=1,bty="n", pch=16, yaxs = "i", xaxs="i", col=colcode)
abline(0,1) 
abline(ctsaph, lty=2)
text(4.5, 8.5, labels="A")

plot(ctlsalog~ photolsalog,xlim=c(4, 8),ylim=c(4,9),
     ylab=expression(paste("Log [CT LSA ", "(mm"^{2},")]")),
     xlab=expression(paste("Log [PH LSA ", "(mm"^{2},")]")), cex=1,yaxs = "i", xaxs="i",bty="n", pch=16, col=colcode)     
abline(0,1) 
graph1<-lm(ctlsalog~ photolsalog)
abline(clsaph, lty=2)
text(4.5, 8.5, labels="B")

plot(ctvollog~ photovollog, xlim=c(2, 9),ylim=c(2,9),
     ylab=expression(paste("Log [CT Vol ", "(mm"^{3},")]")),
     xlab=expression(paste("Log [PH Vol","(mm"^{3},")]")), cex=1,yaxs = "i", xaxs="i",bty="n", pch=16, col=colcode)     
abline(0,1)  
graph1<-lm(ctvollog~ photovollog)
abline(c(cvsaph[1],cvsaph[2]), lty=2, col="red")
abline(c(cvsaph[1]+cvsaph[3],cvsaph[2]+cvsaph[5]), lty=2, col="blue4")
abline(c(cvsaph[1]+cvsaph[4],cvsaph[2]+cvsaph[6]), lty=2, col="chartreuse3")
text(3, 8.5, labels="C")


