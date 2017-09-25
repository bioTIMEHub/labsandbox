# checking what causes discrepancy in results with old vs. new data
idsconstant<-unique(TS$ID)

idsold<-read.csv("final queries/idsold.csv") # list of study ids for previous iteration of analysis
idsoldv<- idsconstant[idsconstant %in% idsold[,2]] # overlap with current study ids 
# (some studies were removed or broken up) note that study 148 in particular is very large and is now absent, check why
idsnewv<- idsconstant[!idsconstant %in% idsold[,2]] # new studies

# make idsconstant equal to each of these and run WandLcode.R 
idsconstant<-idsoldv
idsconstant<-idsnewv

par(mfrow=c(2,2))

hist(as.numeric(colz$slope))
hist(as.numeric(ext$slope))


# these look identical and centered on zero
# repeated with several other rarefactions with same result
# now try with multiple examples of rarefactions from before

TS<- read.csv("final queries/TSrf3.csv")# TSrf has been rarefied for constant number of samples

# extinctions and colonizations have negative and positive slopes 

#cannot see any differences in code for rarefaction
# now testing differences in the raw data 

setwd("~/Dropbox/BIO-TIME/winers-losers/balancedWandL/query")
TS<- read.csv("TSWL.csv")# TS has 178 datasets containing time #series of 10 or more years
IDs <- unique(TS$STUDY_ID)
setwd("~/Dropbox/BIO-TIME/winers-losers old/final queries")

TSold<-read.csv("TS.csv")
IDsold <- unique(TSold$ID)

idscommon<- IDsold[IDsold %in% IDs]

id<-IDsold[1]
nrecords<-c()
tabnd<-c()
for(id in idscommon){
  studynew<- TS[TS$STUDY_ID==id,]
  studyold<- TSold[TSold$ID==id,]
  nrecordsnew<- c(nrecordsnew,dim(studynew)[1])
  nrecordsold<- c(nrecordsold,dim(studyold)[1])
  nrecords<- c(nrecords, dim(studynew)[1]-dim(studyold)[1])
  tabnd<-c(tabnd,sum(studynew$sum.allrawdata.ABUNDANCE)-sum(studyold$Abundance))
}
