#code for running teh assemblage size chaneg model. Within this program i plot the line graph showing slopes
#of change, compare assemblage size chaneg to that of teh null model, abd export assemblage size
#slopes to put into other analysis (teh comparison between assemblage size change and rate of change of dominacne


setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\NullModel_correctsubset")
domAbNull <- read.csv("AssemblageizeNull.csv")

AbundanceDataAll <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\data preperation\\AbudnaceOnceRarefy.csv")
names(AbundanceDataAll)
metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Data_Oct2017\\bioTIMEmetadataJune.csv")

NullStudies <- read.csv("study_IDs.csv")#list of studies used in teh null model 

#compare the mixed model of teh community size (measured in abundance) with the null model results
#that incuded teh years where species where not present 

library(matrixStats)
library(lme4)
library(vegan)
library(paleotree)
library(tidyr)
library(binhf)
library(reshape2)
library(dplyr)
library(ggplot2)

##PREPARING DATA########
----------------------------
#centering each dataset around the mean so that my intercept value will be my mean value as well in the model
AbundanceDataAll$MeanCentredYear <- AbundanceDataAll$Year

for (i in unique(AbundanceDataAll$Study_ID)){
	studyi <- AbundanceDataAll[AbundanceDataAll$Study_ID == i,]	# select data from a singel study
	meani <- mean(studyi$Year)	#calculate mean year
	AbundanceDataAll$MeanYear[AbundanceDataAll$Study_ID == i] <- meani
	meanCentredYeari <- studyi$Year - meani # subtracting the mean from each year 
	AbundanceDataAll$MeanCentredYear [AbundanceDataAll$Study == i] <- meanCentredYeari		#putting data into the dataset
}

head(AbundanceDataAll)

#remove studies with fewer than 10 years
AbSelect <- data.frame(AbundanceDataAll %>%
	group_by(Study_ID)%>%
	mutate(LengthYear = n_distinct(Year))%>%
	filter(LengthYear > 9)
)

length(unique(AbSelect$Study_ID))

#remove biomass studies 
AbStduies <- metaData$STUDY_ID[!is.na(metaData$ABUNDANCE_TYPE)]
rarefiednopaNoBio <- AbSelect[AbSelect$Study_ID %in% AbStduies ,]
length(unique(rarefiednopaNoBio$Study_ID))

#remove presence/absence studies 
metaData[metaData$ABUNDANCE_TYPE == "Presence/Absence",]
paStudies <- metaData[metaData$ABUNDANCE_TYPE == "Presence/Absence", 1]
paStudies <- paStudies[!is.na(paStudies)]
rarefiednopa <- rarefiednopaNoBio[!rarefiednopaNoBio$Study_ID %in% paStudies,]
length(unique(rarefiednopa$Study_ID))

#remove study 360 because although it is listed as having count data, all 
#abundances are 0. Thsi corrected in teh new query
#We also remove study 248 because it is too experimental 

removeSTudies <- c("360", "248")
AbSelect <- rarefiednopa[!rarefiednopa$Study_ID %in% removeSTudies,]

names(AbSelect) 
#running mixed model on actual dominance data 
#--------------------------------------------


#get the assemblage size each year measured in numerical abudnance 
AssSize <- data.frame(AbSelect %>%
		group_by(Study_ID, MeanCentredYear) %>%
		mutate(abSum = sum(Abudance)) %>% #get overall abunacne of all individulas that year/study
		mutate(logAb = log2(1+ abSum)) %>% 
		select(Study_ID, MeanCentredYear, abSum, logAb, MeanYear)
	)

unique(AssSize$Study_ID)
AssSize$Study_ID
head(AssSize )

#run model
AssSizeModel <- lmer( logAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=AssSize)
summary(AssSizeModel)

AssSizeModel_Lines <- coef(AssSizeModel)$Study[2]
AssSizeModel_Lines$Study_ID <- rownames(AssSizeModel_Lines)

#calculate metrics from null model
#-----------------------------------

domAbNull$X <- NULL
domAbNull$Q97.5 <- domAbNull[,1]
domAbNull$Q2.5 <- domAbNull[,1]
names(domAbNull)[1] <- "Study_ID"

names(domAbNull)
for(b in c(1:length(domAbNull$Study_ID))){ #calculates upper and lower quantiles from the null model for each study 1000 repeats 

	q12 <- quantile(domAbNull[b,2:1001],probs=c(.025,.975))[[1]]
	q197 <- quantile(domAbNull[b,2:1001],probs=c(.025,.975))[[2]]
	domAbNull$Q97.5[b] <- q197
	domAbNull$Q2.5[b] <- q12
}

#combine null model with actual resulst 
#------------------------------------------  

domAbNull$Study_ID <- as.character(domAbNull$Study_ID)
AssSizeModel_Lines$Study_ID <- as.character(AssSizeModel_Lines$Study_ID)

AssNullAndActual<- merge(domAbNull, AssSizeModel_Lines, by = "Study_ID")#incuding the slopes values from the actual mixed model fo data
names(AssNullAndActual)[1004] <- "ActualSlopes"

AssNullAndActual[c(1:5),c(1000:1004)]

#how many studies show a noticable change?
#---------------------------------------------

#selecting studies where the dominant abundance changes more than woudl be expected from the null model
AssNullAndActual$AssemblageSizeChange<-AssNullAndActual$Study_ID #making a spare column to input data 

AssNullAndActual$AssemblageSizeChange[AssNullAndActual$ActualSlopes < AssNullAndActual$Q2.5] <- "decrease"
AssNullAndActual$AssemblageSizeChange[AssNullAndActual$ActualSlopes > AssNullAndActual$Q97.5] <- "increase"
AssNullAndActual$AssemblageSizeChange[AssNullAndActual$ActualSlopes < AssNullAndActual$Q97.5 & AssNullAndActual$ActualSlopes > AssNullAndActual$Q2.5] <- "noChange" 

barplot(table(AssNullAndActual$AssemblageSizeChange))


#get actual percentiles for plotting
#----------------------------------------

percentilleData <- data.frame(matrix(NA, length(domAbNull$Study_ID[!domAbNull$Study_ID %in% removeSTudies]), 2))
names(percentilleData) <- c("Study_ID", "AssSize")
percentilleData[,1] <- domAbNull$Study_ID[!domAbNull$Study_ID %in% removeSTudies]

for (p in c(1:length(domAbNull$Study_ID))){

	pCurve <- ecdf(domAbNull[p,2:1001])(AssNullAndActual$ActualSlopes[p])
	percentilleData[p, 2] <- pCurve
}

names(AssNullAndActual)
merge(percentilleData, AssNullAndActual[,c(1, 1005)], by = "Study_ID")

#write.csv(percentilleData, "AssSizeChange_Percentiles.csv")


#calculate Z scores : (slope - meanslope)/standard deviation slope 
###-------------------------

AssNullAndActual$nullmean <- rowMeans(AssNullAndActual[,2:1001])
names(AssNullAndActual)
AssNullAndActual$SD <- apply(AssNullAndActual[,2:1001],1,sd)

AssNullAndActual$AssSizeZscores <- (AssNullAndActual$ActualSlopes - AssNullAndActual$nullmean)/AssNullAndActual$SD

#not working yet
#write.csv(AssNullAndActual[,c(1, 1008)], "AssSize_Zscores.csv")

###making a sloep plot for presentations
#----------------------------------------------
names(AssSize)
head(metaData)
aDomLines <- coef(AssSizeModel)$Study_ID
aDomLines$Study <- rownames(aDomLines) #extract the ranom effect levels
lines2 <- merge(aDomLines, metaData[,c(1, 15, 14)], by.x = "Study", by.y = "STUDY_ID", all = FALSE) #add this infor to info on line slopes 
lines <- merge(lines2, AbSelect[,c(5,8,6)], by.x = "Study", by.y = "Study_ID")


lines <-lines[!duplicated(lines), ]
names(lines)[2] <- "intercept"
names(lines)[3] <- "slope"
head(lines)

#work out start and end point of my studies (random effects)  
lines$meanCentedStartYear <- lines$START_YEAR - lines$MeanYear
lines$meanCentedEndYear <- lines$END_YEAR - lines$MeanYear
lines$StartYActuall <- lines$slope * lines$meanCentedStartYear + lines$intercept
lines$EndYActuall <- lines$slope * lines$meanCentedEndYear + lines$intercept
AssSize$ActualYear <- AssSize$MeanYear + AssSize$MeanCentredYear

#work out start and end point of my studies of actuall dates (random effects)  
lines$StartY <- lines$slope * lines$meanCentedStartYear + lines$intercept
lines$EndY <- lines$slope * lines$meanCentedEndYear + lines$intercept

summary(DominancModel)
#manualy calculate the mean slope
endypoint <- 0.0069 * 127 + 13.4

head(domAb)

#plot the same but with points see through and no colour differences 
m3.plot2 <- ggplot(data = AssSize, aes(x = ActualYear, y = logAb)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plotA2 <- m3.plot2 + theme_classic()+
	theme(text = element_text(size=30)) +
	geom_segment(data =  lines, aes(x = START_YEAR, y = lines$StartY, xend = lines$END_YEAR, yend = lines$EndY),alpha = 2/3, colour = "#333333")+
	geom_segment(aes(x = 1899, y = 12.36, xend = 2017, yend = endypoint), colour = "blue",lwd=1.25) + 
	labs(x="Year",y="Log Assemblge Size")+
	scale_x_continuous(limits = c(1899, 2017))+ 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank()) 


