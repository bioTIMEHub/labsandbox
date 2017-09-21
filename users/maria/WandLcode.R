TS<- read.csv("TSrfFiles/TSrfwinlLosev1.csv")# TSrf has been rarefied for constant number of samples

nsampless<-list()
idplace<-1
for(id in unique(TS$ID)){
  Year<-TS$Year[TS$ID==id]
  SampleID<-TS$SampleID[TS$ID==id]
  nsamples<-c()
  for(y in unique(Year)){
    nsamples<-c(nsamples, length(unique(SampleID[Year==y])))
  }
  nsampless[[idplace]]<-nsamples
  idplace<-idplace+1
}
ids<-unique(TS$ID)
names(nsampless)<-ids
#find time series with constant number of samples
maxes<-mapply(FUN=max, nsampless)
mins<-mapply(FUN=min, nsampless)
constant<-maxes-mins

idsconstant<-ids[constant==0]# checking that all ids have constant number of samples


#####################################################################################
# functions


library(tseries) #make sure the tseries package is installed
get.coeff <- function(x,y) coef(summary(lm(y~x)))[2,c(1,4)] # function for slope and p

sig.tests <- function(vec=runif(100)){
  
  vec <- sort(vec)  
  
  # simple significance, a count of tests with p <0.05 
  simple <- length(vec[vec <= 0.05])
  
  # classic bonferroni with 0.05 level divided by number of tests
  bonfer <- length(vec[vec <= 0.05/length(vec)]) 
  
  # bh Benjamini and Hochberg (1995) test for independent or 
  # positively correlated families of tests
  #--------------------------------------------------------
  # Benjamini, Y. and Y. Hochberg. 1995. Controlling the 
  # False Discovery Rate: a practical and powerful approach
  # to multiple testing.
  # J. Royal Stat Soc Ser B 57: 289-300.
  #-------------------------------------------------------
  comp <- 0.05*(1:length(vec))/length(vec)
  bh <- sum(vec <= comp)
  
  # bl Benjamini and Liu (1999) test for positively and negatively
  # covarying tests.
  
  # NOTE: With randomly constructed toy data sets in which all the tests
  # are independent, the bl test often comes out too
  # conservative, and seems to match more closely the classic bonferroni correction.
  # The bh performs much better and gives surprisingly accurate results for 
  # data sets that consist of a mixture of significant and non-significant tests.
  
  # For the winners and losers analysis, it might be best to use the more liberal
  # bh test, since we want to be able to detect things that are there.
  
  # It would be best in the winners losers analysis to apply the test separately
  # to each time series, because this will give the best chance of finding something
  # in each individual series. It probably will be too conservative to lump data
  # across studies.
  # ------------------------------------------------------
  # Benjamini, Y. and W. Liu. 1999. A distribution-free multiple test procedure 
  # that controls the false discovery rate.Tel Tel Aviv, RP-SOR-99-3:
  # Department of Statistic and O.R. Tel Aviv University.
  
  # Could also cite instead
  # Benjamini, Y., D. Drai, G. Elmer, N. Kafkafi, and I. Golani. 2001.
  # Controlling the false discovery rate in behavior genetics research.
  # Behavioural Brain Research 125: 279-284.
  # NOTE: Fully-worked example in Table 1 of this paper was used to confirm
  # correct output of this function.
  # ------------------------------------------------------
  comp1 <- pmin(rep(0.05,length(vec)),((0.05*length(vec))/((length(vec) + 1 - (1:length(vec))))^2))
  
  bl <- sum(vec <= comp1)
  return(list(simple=simple,bonfer=bonfer,bh=bh,bl=bl,raw=vec[1:simple]))
}

N.turnovers <- function (vec=rbinom(50,1,0.5)) {
  
  z <- diff(vec)   # get difference lag one
  col <- sum(z==1) # count 0 -> 1
  ext <- sum(z==-1)# count 1 -> 0
  
  result <- c(col,ext)
  names(result) <- c("col","ext")
  return(result)
  
}


#########################################################################

#run winer and loser analysis for studies with constant sampling effort

WLEC<-data.frame(ID=0,Species=0,runs.pvalue=0,col=0,ext=0,slope=0,slope.pvalue=0,Ninit=0,meanN=0)
idplace<-1
for(id in idsconstant){
  # getting data for relevant studyID
  data<-TS[TS$ID==id,]
  data<- data[data$Abundance>0,]
  groups<-data.frame(as.character(data$Species),as.numeric(data$Year))
  data.mat<- tapply(data$Abundance,groups, FUN=sum)
  # formating data into species by time matrix
  data.mat[is.na(data.mat)]<-0
  #removing species that are always absent
  N.species<-dim(data.mat)[1] #getting number of species
  WLEC[idplace:(idplace+N.species-1),]<-cbind(rep(id,N.species),rownames(data.mat),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),rep(NA,N.species),data.mat[,1],apply(data.mat,1,mean))
  
  ############### identifying extinctions and colonizations
  Binary.Data <- data.mat #using previously simulated data matrix
  Binary.Data[Binary.Data > 0] <- 1 # convert to binary
  
  # remove rows in which there were no absences
  Binary.Data <- Binary.Data[which(rowSums(Binary.Data)< ncol(Binary.Data)),]
  
  # create data frame to hold results
  if(!is.matrix(Binary.Data)) Binary.Data <- t(Binary.Data)# HS added
  
  if(dim(Binary.Data)[1]>0){
    # loop through the data
    for (i in 1:nrow(Binary.Data)) {
      # Extract data for a species, conduct runs test, save output
      z <- Binary.Data[i,]
      runs.p <- runs.test(as.factor(z),alternative="less")
      # Count number of colonizations and extinctions, save output
      turnover <- N.turnovers(z)
      WLEC[WLEC$ID==id & WLEC$Species==rownames(Binary.Data)[i],3:5]<-c(runs.p$p.value, turnover[1], turnover[2])
    }
  }
  ##### Winers and Losers
  data.mat<-t(apply(data.mat,1,scale))# scaling by subtracting mean and dividing by sd
  Time<- unique(data$Year) #getting time vector  
  ######### identifying species with positive (Winers) and negative (Losers) trends
  WLEC[idplace:(idplace+N.species-1),6:7]<-t(apply(data.mat,1,get.coeff,x=Time)) # apply to each row of data

  idplace<-idplace+N.species
}




################################################
# replacing slopes of time series where an extinction or colonization occurred
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


placemixed<-which(WLEC$runs.pvalue<0.05)
placemixed<-subset(placemixed,!(placemixed %in%c(placeext,placecol)))
WLEC$runs.pvalue<-as.numeric(WLEC$runs.pvalue)
WLEC$slope<-as.numeric(WLEC$slope)
WLEC$slope.pvalue<-as.numeric(WLEC$slope.pvalue)
write.csv(WLEC,"WLECFiles/WLECfinalv1.csv")
#############################################
nozeros<-subset(WLEC, is.na(runs.pvalue))
ext<- subset(WLEC, col==0&ext==1 & runs.pvalue<0.05)
colz<-subset(WLEC, col==1&ext==0 & runs.pvalue<0.05)
nr0<-WLEC[placemixed,]
r0<-subset(WLEC, runs.pvalue>0.05)


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
