# checking what causes discrepancy in results with old vs. new data
idsconstant<-unique(TS$ID)

idsold<-read.csv("final queries/idsold.csv") # list of study ids for previous iteration of analysis
idsoldv<- idsconstant[idsconstant %in% idsold[,2]] # overlap with current study ids (some studies were removed or broken up)
idsnewv<- idsconstant[!idsconstant %in% idsold[,2]] # new studies

# make idsconstan! equal to each of these and run WandLcode.R 
idsconstant<-idsoldv
idsconstant<-idsnewv
# these look identical and centered on zero
# repeated with several other rarefactions with same result
par(mfrow=c(2,2))

hist(as.numeric(colz$slope))
hist(as.numeric(ext$slope))

idsconstant<-idsnewv


TS<- read.csv("final queries/TSrf3.csv")# TSrf has been rarefied for constant number of samples
