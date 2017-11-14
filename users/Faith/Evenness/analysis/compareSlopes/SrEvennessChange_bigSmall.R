setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
library(ggplot2)
library(ggExtra)
library(dplyr)

###thsi code is to see if there is a different corrilation between evenness change 
and species rcihness change in small assemblages (>20) vs large assemblages 

evennessSim <- read.csv("evennessChangeSimp.csv")
evenness <- read.csv("evennessChange.csv")
SRChange <- read.csv("SpeciesRichnessChange.csv")
RarefiedData <- read.csv("AbudnaceOnceRarefy.csv")
head(RarefiedData)

##preparying data
names(evennessSim)[4] <- "EvennessSim1"
names(evenness)[4] <- "Evenness"

##merging all data and selecting columns of interest
SrEv <- merge(SRChange, evenness, by.x = "Study_ID", by.y = "Study")
AllData <- merge(SrEv, evennessSim, by.x = "Study_ID", by.y = "Study")

names(AllData)
SelectData <- AllData[,c(1,4,9,10,18)]
head(SelectData)

##calculate mean and mimimum species richness numbers to seperate
##out assemblages with low species richness

#mimimum species richness
SrMin<- RarefiedData %>% 
	group_by(Study_ID, Year) %>%
	summarise(SRNumber = n())%>%
	summarise(SRMin = min(SRNumber))
MinimumSr <- data.frame(SrMin)
minSmall <- MinimumSr[MinimumSr$SRMin < 20,]

#mean species richness 
SrMean<- RarefiedData %>% 
	group_by(Study_ID, Year) %>%
	summarise(SRNumber = n())%>%
	summarise(SRMean = mean(SRNumber))

MeanSr <- data.frame(SrMean)
meanSmall <- MeanSr[MeanSr$SRMean < 20,]

#merge slopes with minimum species richness data
SelecDataMinSize <- merge(SelectData, MinimumSr , by = "Study_ID")


##select data with assemblages species richness of less than 20 
#to look at correlation 
smallSR <- SelecDataMinSize[SelecDataMinSize$SRMin < 20, ]
BigSR <- SelecDataMinSize[SelecDataMinSize$SRMin > 20, ]
cor(smallSR$SRslope, smallSR$EvennessSim1 ,method = "spearman")
cor(BigSR$SRslope, BigSR$EvennessSim1 ,method = "spearman")

#Plotting results 
#------------------------
bigSR <- ggplot(BigSR, aes(EvennessSim1, SRslope)) + 
	geom_point() + 
	theme_classic()+ 
	labs(x = "Evenness Change", y = "Species Richness Change") +
	scale_x_continuous(limits = c(-0.009,0.009))+
	scale_y_continuous(limits = c(-0.09,0.09))
ggMarginal(bigSR , type = "histogram")

LittleSR <- ggplot(smallSR, aes(EvennessSim1, SRslope)) + 
	geom_point() + 
	theme_classic()+ 
	labs(x = "Evenness Change", y = "Species Richness Change")+
	scale_x_continuous(limits = c(-0.009,0.009))+
	scale_y_continuous(limits = c(-0.09,0.09))
ggMarginal(LittleSR , type = "histogram")









