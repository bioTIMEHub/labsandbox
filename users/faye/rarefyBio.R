TS<- read.csv("TS.csv")# TS has 178 datasets containing time #series of 10 or more years


### from Maria's code

rarefysamples<-function(Year, SampleID, Species, Abundance, resamps) {
  #######################################################################
  # takes as input a  Year, SampleID, Species, Abundance and number of resamples
  # which should be in dataframe so that elements match
  # calculates turnover:
  # 1) between each year and the first year 
  # 2) between pairs of adjacent years 
  # 3) between each year and the lasy year of the time series
  # for the rarefied pooled samples
  ###########################################################################
  
  rareftab<-data.frame(array(NA,dim=c(0,3)))
  # getting vector with number of samples per year
  nsamples<-c()
  for(y in unique(Year)){
    nsamples<-c(nsamples, length(unique(SampleID[Year==y])))
  }
  t<-1
  minsample<-min(nsamples)
  for(repeats in 1:resamps){
    raref<-data.frame(array(NA,dim=c(1,3)))
    names(raref)<-c("Year","Species","Abundance")
    for(y in unique(Year)){
      #getting samples for this year
      samps<-unique(SampleID[Year==y])
      # re-sampling samples to equalize number of samples
      sam<-as.character(sample(samps,minsample,replace=T))
      # getting data that belongs to bootstraped samples
      rarefyear<-data.frame(SampleID[which(SampleID %in% sam & Year == y)],Species[which(SampleID %in% sam & Year == y)],Abundance[which(SampleID %in% sam & Year == y)])
      names(rarefyear)<-c("SampleID", "Species", "Abundance")
      # calculating pooled abundances of each species to store
      spabun<-tapply(as.numeric(rarefyear[,3]),as.character(rarefyear[,2]),sum)
      spar<-data.frame(rep(y, length(spabun)),names(spabun),spabun, row.names=NULL)
      names(spar)<-c("Year","Species","Abundance")
      raref<-rbind(raref,spar)
    }
    # calculating year by species table of abundance
    rareftab<-rbind(rareftab,cbind(rep(repeats,dim(raref)[1]),raref))
  }
  return(rareftab)
}


####################end of function###############


TSrf <- list()
IDs <- unique(TS$STUDY_ID)
#tu use function
for(i in 1:length(IDs)){
data<-TS[TS$STUDY_ID==IDs[i],]
TSrf[[i]]<-rarefysamples(data$YEAR, data$SAMPLE_DESC, data$SPECIES, data$ABUNDANCE, 1)
}
names(TSrf) <- IDs

test <- do.call(rbind, TSrf)
test <- data.frame(test, ID=rep(names(TSrf), times=unlist(lapply(TSrf, nrow))))
test <- test[!is.na(test$YEAR),-1]

write.csv(test, file="TsrfwinLoseV2.csv", row.names=F) #change name appropriately

