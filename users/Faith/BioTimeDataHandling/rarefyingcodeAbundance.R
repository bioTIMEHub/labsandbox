setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
library(lmtest)
library(plyr)
library(ggplot2)
library(rms)
library(nlme)
library(mgcv)
library(segmented)
library(grid)
library(gridExtra)
library(scales)
library(sme)
library(parallel)
library(foreach)
library(doParallel)

# Maria Dornelas 13.05.2013
# rarefying for equal number of samples

rarefysamplesSNE<-function(Year, SampleID, Species, Abundance, resamps) {
  #######################################################################
  # takes as input a  Year, SampleID and and Species 
  # which should be in dataframe so that elements match
  # must have Chao estimator function on workspace
  ###########################################################################
  library(vegan)
  
  ts<-data.frame(array(NA,dim=c(length(unique(Year)),13)))
  names(ts)<-c("year", "S","S2","varS","N","N2", "varN", "SsqrtN","PIE","DomMc","expShannon", "Chao", "Chao2")
  # getting vector with number of samples per year
  nsamples<-c()
  for(y in unique(Year)){
    nsamples<-c(nsamples, length(unique(SampleID[Year==y])))
  }
  t<-1
  minsample<-min(nsamples)
  for(y in unique(Year)){
    rarefS<-c()
    rarefS2<-c()
    rarefvarS<-c()
    rarefN<-c()
    rarefN2<-c()
    rarefvarN<-c()
    rarefSN<-c()
    rarefPIE<-c()
    rarefDomMc<-c()
    rarefexpShan<-c()
    rarefChao<-c()
    rarefChao2<-c()
    #getting samples for this year
    samps<-unique(SampleID[Year==y])
    for(repeats in 1:resamps){
      # re-sampling samples to equalize number of samples
      sam<-as.character(sample(samps,minsample, replace=T))
      raref<-cbind(SampleID[which(SampleID %in% sam & Year == y)],Species[which(SampleID %in% sam & Year == y)],Abundance[which(SampleID %in% sam & Year == y)])
      # calculating pooled abundance of each species in re-sampled samples
      sad<-tapply(as.numeric(raref[,3]),raref[,2],sum)
      sad<-sort(sad, decreasing=T)
      # number of species per sample
      Spersample<-tapply(raref[,2],raref[,1],length)
      # total abundance per sample
      Npersample<-tapply(as.numeric(raref[,3]),raref[,1],sum)
      # calculating number of species in re-sampled samples
      rarefS<-c(rarefS, length(unique(Species[which(SampleID %in% sam & Year == y)])))
      rarefS2<-c(rarefS2, length(sad))
      rarefvarS<-c(rarefvarS,var(Spersample))
      rarefN<-c(rarefN, sum(Abundance[which(SampleID %in% sam & Year == y)]))
      rarefN2<-c(rarefN2, sum(sad))
      rarefvarN<-c(rarefvarN,var(Npersample))
      rarefSN<-c(rarefSN, length(sad)/(sum(sad)^(1/2)))
      rarefPIE<-c(rarefPIE, (sum(sad)/(sum(sad)-1))*(1-(sum((sad/sum(sad))^2))))#Hurlbert 1971
      rarefDomMc<-c(rarefDomMc,as.vector((sad[1]+sad[2])/sum(sad)))
      rarefexpShan<-c(rarefexpShan,exp(diversity(sad,index="shannon")))
      rarefChao<-c(rarefChao,Chao(as.vector(sad)))
      rarefChao2<-c(rarefChao2,Chao2(as.vector(sad)))
    }
    ts[t,]<-c(y,median(rarefS, na.rm=T),median(rarefS2, na.rm=T),median(rarefvarS, na.rm=T),median(rarefN, na.rm=T),median(rarefN2, na.rm=T),median(rarefvarN, na.rm=T), median(rarefSN, na.rm=T),median(rarefPIE, na.rm=T),median(rarefDomMc, na.rm=T),median(rarefexpShan, na.rm=T), median(rarefChao, na.rm=T), median(rarefChao2, na.rm=T))
    t<-t+1
  }
  return(ts)
}

########################################
Chao<-function(sad){
  nsing<-sum(sad==1)
  ndoub<-sum(sad==2)
  if(ndoub==0){ Chao<-NA}
  else{ 
    S<-length(sad)
    Chao<-S+nsing^2/ndoub/2}
  return(Chao)
}

########################################
# Chao1 function
# takes as input a vector of species abundance counts (including 0)
# returns bias-corrected Chao1 (Eq. 2 in EstimateS manual)
Chao2 <- function(v=rpois(20,1))
{
  Sobs <- length(v[v>0])
  F1 <- length(v[v==1])
  F2 <- length(v[v==2])
  Chao2 <- Sobs +(F1*(F1-1))/(2*(F2+1))
  return(Chao2)
}

#
#--------------------------------------------------------------------------------------------------end of maria's code
#

#-------------------------------------------
#Reading in Abundance data and rarefying data
#-------------------------------------------

DataNoSmall <- read.csv("AbundanceDataNoSmallVales.csv")
head(DataNoSmall)
#asking r to use 3 cores

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)

cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)

#rarefying

DataRare <- foreach (s = unique(DataNoSmall$STUDY_ID)) %dopar% {
	datai <- DataNoSmall[DataNoSmall$STUDY_ID == s,]
	rareData <- rarefysamplesSNE(Year = datai$YEAR, SampleID = datai$SAMPLE_DESC, Species = datai$GENUS_SPECIES, Abundance = datai$ABUNDANCE, resamps= 200)
	rareData$Study <- s
	rareData

}

stopCluster(cl)

newDFAll <- do.call(rbind, DataRare)#bind all lists together
row.names(newDFAll) <- NULL 
head(newDFAll)
write.csv(newDFAll, "rarefiedAbundanceNoSmall.csv")


