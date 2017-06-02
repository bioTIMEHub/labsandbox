# population simulations to compare histograms of slopes and geometric mean index
# Maria Dornelas 2 June 2017
# Work in progress
# functions ---------------------------------------------------------------
# function to simulate exponential populations
expop<- function(r,Ninit,tl){
  # exponential population growth model
  # r is the instantaneous rate of increase, Ninit is the initial population size
  # tl is the number of time steps in population simulation
  Ns<-Ninit*exp(r*1:tl)
  Ns2<-Ninit
  return(Ns)
}

# Nick's functiont o get slope and p value of OLS linear trend of pop size and time
get.coeff <- function(x,y) coef(summary(lm(y~x)))[2,c(1,4)] # function for slope and p

# simulations -------------------------------------------------------------
# simulations for 100 pops with positive r
rs<-runif(100,0,0.03)
Ninits<-runif(100,1,100)
Np<-matrix(NA,100,50)
for(pops in 1:100){
  Np[pops,]<-expop(rs[pops],Ninits[pops],50)
}  

# getting slopes for each population

Poptrends<-t(apply(Np,1,get.coeff,x=1:50)) # apply to each row of data
hist(Poptrends[,1])


# LPI ---------------------------------------------------------------------

# this bit is not working
for(pops in 1:100){
LogN[[pops]] <- data.frame(Year=1:50, log(Np[pops]))
} 

d.bar <- lapply(LogN, function(x)apply(apply(x[,-1], 2, diff), 1, mean))
LPI <- lapply(d.bar, function(x)c(1, 10^cumsum(x)))

LPI2 <- list()
for(i in 1:length(LPI)){
  LPI2[[i]] <- data.frame(Year=LogN[[i]][,1], dbar=LPI[[i]])
}

