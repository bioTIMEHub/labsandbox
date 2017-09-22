# Prepared by Maria Dornelas 21 september 2017, using some functions prepared by Nick Gotelli

# this script takes an object with sample rarefied species abundances for bioTIME studies with more
# than 10 years in duration. 
# For each population within each study, we:
# - identify colonizations and extinctions using a runs test
# - estimate linear slopes 
# - estimate linear slopes that use only time series before extinctions or after colonizations

# cleaned up older code prepared for the community regulation Science Advances paper
# includes only code relevant for the winners and loosers analysis

rm(list=ls())

TS<- read.csv("TSrfFiles/TSrfwinlLosev1.csv")# TSrf has been rarefied for constant number of samples

idsconstant<-unique(TS$ID)

# functions #########################################################################

library(tseries) #make sure the tseries package is installed
get.coeff <- function(x,y) coef(summary(lm(y~x)))[2,c(1,4)] # function for slope and p

N.turnovers <- function (vec=rbinom(50,1,0.5)) {
  
  z <- diff(vec)   # get difference lag one
  col <- sum(z==1) # count 0 -> 1
  ext <- sum(z==-1)# count 1 -> 0
  
  result <- c(col,ext)
  names(result) <- c("col","ext")
  return(result)
  
}


## loop across pops#############################################################

#run winner and loser analysis

WLEC<-data.frame(ID=0,Species=0,runs.pvalue=0,col=0,ext=0,slope=0,slope.pvalue=0,Ninit=0,meanN=0)
idplace<-1
for(id in idsconstant){
  # getting data for relevant studyID 
  data<-TS[TS$ID==id,]
  groups<-data.frame(as.character(data$Species),as.numeric(data$Year))
  data.mat<- tapply(data$Abundance,groups, FUN=sum)
  # formating data into species by time matrix 
  data.mat[is.na(data.mat)]<-0
  N.species<-dim(data.mat)[1] #getting number of species
  WLEC[idplace:(idplace+N.species-1),]<-cbind(rep(id,N.species),rownames(data.mat),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),data.mat[,1],apply(data.mat,1,mean))
  
  # identifying extinctions and colonizations ####
  Binary.Data <- data.mat #using previously simulated data matrix
  Binary.Data[Binary.Data > 0] <- 1 # convert to binary
  
  # remove rows in which there were no absences
  Binary.Data <- Binary.Data[which(rowSums(Binary.Data)< ncol(Binary.Data)),]
  
  # create data frame to hold results
  if(!is.matrix(Binary.Data)) Binary.Data <- t(Binary.Data)# HS added
  
  if(dim(Binary.Data)[1]>0){
    # loop runs test through the data
    for (i in 1:nrow(Binary.Data)) {
      # Extract data for a species, conduct runs test, save output
      z <- Binary.Data[i,]
      runs.p <- runs.test(as.factor(z),alternative="less")
      # Count number of colonizations and extinctions, save output
      turnover <- N.turnovers(z)
      WLEC[WLEC$ID==id & WLEC$Species==rownames(Binary.Data)[i],3:5]<-c(runs.p$p.value, turnover[1], turnover[2])
    }
  }

  # identifying Winners and Loosers ########################################################
  data.mat<-t(apply(data.mat,1,scale))# scaling by subtracting mean and dividing by sd
  Time<- unique(data$Year) #getting time vector  
  # identifying species with positive (Winners) and negative (Losers) trends
  WLEC[idplace:(idplace+N.species-1),6:7]<-t(apply(data.mat,1,get.coeff,x=Time)) # apply to each row of data

  idplace<-idplace+N.species
}

# replacing slopes ################################################
# of time series where an extinction or colonization occurred
# with slopes estimated for when the species was present only

ltsext<-c()
placeext<-which(WLEC$col==0 & WLEC$ext==1 & WLEC$runs.pvalue<0.05)
pl<-1
for(i in placeext){
  data<-TS[TS$ID==WLEC$ID[i] & TS$Species==WLEC$Species[i],]
  ltsext<-c(ltsext,dim(data)[1])
  if(ltsext[pl]>3){
    WLEC[i,6:7]<-coef(summary(lm(data$Abundance~data$Year)))[2,c(1,4)]
    pl<-pl+1
  }
}
  

ltscol<-c()
placecol<-which(WLEC$col==1 & WLEC$ext==0 & WLEC$runs.pvalue<0.05)
pl<-1
for(i in placecol){
  data<-TS[TS$ID==WLEC$ID[i] & TS$Species==WLEC$Species[i],]
  ltscol<-c(ltscol,dim(data)[1])
  if(ltscol[pl]>3){
    WLEC[i,6:7]<-coef(summary(lm(data$Abundance~data$Year)))[2,c(1,4)]
    pl<-pl+1
  }
}

# saving results ####

placemixed<-which(WLEC$runs.pvalue<0.05)
placemixed<-subset(placemixed,!(placemixed %in%c(placeext,placecol)))
WLEC$runs.pvalue<-as.numeric(WLEC$runs.pvalue)
WLEC$slope<-as.numeric(WLEC$slope)
WLEC$slope.pvalue<-as.numeric(WLEC$slope.pvalue)
write.csv(WLEC,"WLECFiles/WLECfinalv1.csv")

###subsetting results###################################

nozeros<-subset(WLEC, is.na(runs.pvalue)) #trends for populations always present 
ext<- subset(WLEC, col==0&ext==1 & runs.pvalue<0.05) # trends for pops that go extinct
colz<-subset(WLEC, col==1&ext==0 & runs.pvalue<0.05) # trends for new colonizations
nr0<-WLEC[placemixed,] # trends for non random multiple extinctions and colonizations
r0<-subset(WLEC, runs.pvalue>0.05) # trends for random extinctions and colonizations


negslope<- subset(WLEC, slope.pvalue<0.05 & slope<0)
posslope<- subset(WLEC, slope.pvalue<0.05 & slope>0)
zslope<- subset(WLEC, slope.pvalue>0.05)

negslopeext<- subset(ext, slope.pvalue<0.05 & slope<0)
posslopeext<- subset(ext, slope.pvalue<0.05 & slope>0)
zslopeext<- subset(ext, slope.pvalue>0.05)

negslopecol<- subset(colz,slope.pvalue<0.05 & slope<0)
posslopecol<- subset(colz,slope.pvalue<0.05 & slope>0)
zslopecol<- subset(colz,slope.pvalue>0.05)

negslopenoz<- subset(nozeros,slope.pvalue<0.05 & slope<0)
posslopenoz<- subset(nozeros,slope.pvalue<0.05 & slope>0)
zslopenoz<- subset(nozeros,slope.pvalue>0.05)

negsloper0<- subset(r0,slope.pvalue<0.05 & slope<0)
possloper0<- subset(r0,slope.pvalue<0.05 & slope>0)
zsloper0<- subset(r0,slope.pvalue>0.05)

negslopenr0<- subset(nr0,slope.pvalue<0.05 & slope<0)
posslopenr0<- subset(nr0,slope.pvalue<0.05 & slope>0)
zslopenr0<- subset(nr0,slope.pvalue>0.05)

par(mfrow=c(2,3))
hist(as.numeric(nozeros$slope))
hist(as.numeric(r0$slope))
hist(as.numeric(nr0$slope))
hist(as.numeric(colz$slope))
hist(as.numeric(ext$slope))
hist(as.numeric(WLEC$slope))
