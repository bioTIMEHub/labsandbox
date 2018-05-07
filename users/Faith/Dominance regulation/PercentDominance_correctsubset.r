#This code anylises change in relative (%) dominance. various mixed model sare run to see which model fits the best
# It also looks at the empirical data vs the null model results  


setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\NullModel_correctsubset")
domAbNull <- read.csv("DominantPercentNull.csv")
domAbNull[,1]

NullStudies <- read.csv("study_IDs.csv")#list of studies used in teh null model 

AbundanceDataAll <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\data preperation\\AbudnaceOnceRarefy.csv")
names(AbundanceDataAll)
metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Data_Oct2017\\bioTIMEmetadataJune.csv")

#compare the mixed model of relative (%) abundance of dominant species with the null model results
#that incuded teh years hwere species where not present 

library(matrixStats)
library(lme4)
library(vegan)
library(paleotree)
library(tidyr)
library(binhf)
library(reshape2)
library(dplyr)
library(ggplot2)


AbundanceDataAll$Study_ID<- as.character(AbundanceDataAll$Study_ID)


###this code compares the output of the null model against the actual data results for the abundance of dominant species


head(AbundanceDataAll) 

#centering each dataset around the mean so that my intercept value will be my mean value as well in the model
AbundanceDataAll$MeanCentredYear <- AbundanceDataAll$Year

#get mean centred years by considering years where data was not sampled too 
AbundanceDataAll2 <- data.frame(AbundanceDataAll %>% 
	group_by(Study_ID)%>%
	mutate(StartYear = min(Year), endYear = max(Year), nyear= endYear - StartYear) %>%
	mutate(MeanYear = StartYear + nyear/2, MeanCentredYear = Year - MeanYear))

head(AbundanceDataAll2)

AbundanceDataAll2$MeanYear <- AbundanceDataAll2$Abudance

for(study in unique(AbundanceDataAll2$Study_ID)) {
	studyDatai <- AbundanceDataAll2[AbundanceDataAll2$Study_ID == study,]
	allYears <- studyDatai[1,7]:studyDatai[1,8]
	meanYear <- mean(allYears)
	AbundanceDataAll2$MeanYear[AbundanceDataAll2$Study_ID == study] <- meanYear
}


#count length of studies
AbSelect <- data.frame(AbundanceDataAll2 %>%
	group_by(Study_ID)%>%
	mutate(LengthYear = n_distinct(Year))%>%
	filter(LengthYear > 9)
)

#remove biomass studies 
AbStduies <- metaData$STUDY_ID[!is.na(metaData$ABUNDANCE_TYPE)]
rarefiednopaNoBio <- AbSelect[AbSelect$Study_ID %in% AbStduies ,]

#remove presence/absence studies 
metaData[metaData$ABUNDANCE_TYPE == "Presence/Absence",]
paStudies <- metaData[metaData$ABUNDANCE_TYPE == "Presence/Absence", 1]
paStudies <- paStudies[!is.na(paStudies)]
rarefiednopaS <- rarefiednopaNoBio[!rarefiednopaNoBio$Study_ID %in% paStudies,]

#remove study 360 because although it is listed as having count data, all 
#abundances are 0. Thsi corrected in teh new query
#We also remove study 248 because it is too experimental 

removeSTudies <- c("360", "248")
rarefiednopaS <- rarefiednopa[!rarefiednopa$Study_ID %in% removeSTudies,]

head(rarefiednopa)

#get the dominant species each year and output its abundance and percentage abundance
domAbS <- data.frame(rarefiednopaS %>%
		group_by(Study_ID, MeanCentredYear) %>%
		mutate(abSum = sum(Abudance)) %>% #get overall abunacne of all individulas that year/study
		filter(Abudance == max(Abudance)[1]) %>% #select the largest value (if there is more than one tied then just chose the first value)
		filter(1:n() == 1) %>% #make sure  there is only oen value!
		mutate(logDomAb = log2(Abudance + 1), #calculate log abundance
			percentAb = (Abudance/abSum * 100)) #calculate percentage abundance
	)

head(domAbS)
length(unique(domAb$Study_ID))

#remove study 360 because although it is listed as having count data, all 
#abundances are 0. Thsi corrected in teh new query
domAb <- domAbS[!domAbS$Study_ID %in% removeSTudies,]

#run model

DominancModel <- lmer(percentAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAb)
summary(DominancModel)

DominancModel_Lines <- coef(DominancModel)$Study[2]
DominancModel_Lines$Study_ID <- rownames(DominancModel_Lines)
#calculate metrics from null model
#-----------------------------------

domAbNull$X 
NullStudies$number <- rownames(NullStudies)
domAbNullA <- merge(domAbNull, NullStudies, by.x = "X", by.y = "number")

#remove problem studies
domAbNullC <- domAbNullA[!domAbNullA$V1 == "248",]

domAbNullC$Q97.5 <- domAbNullC[,1]
domAbNullC$Q2.5 <- domAbNullC[,1]
names(domAbNullC)[1003] <- "Study_ID"

domAbNullC$Study_ID

for(b in c(1:length(domAbNullC$Study_ID))){ #calculates upper and lower quantiles from the null model for each study 1000 repeats 

	q12 <- quantile(domAbNullC[b,2:1001],probs=c(.025,.975))[[1]]
	q197 <- quantile(domAbNullC[b,2:1001],probs=c(.025,.975))[[2]]
	domAbNullC$Q97.5[b] <- q197
	domAbNullC$Q2.5[b] <- q12
}

#combine null model with actual resulst 
#------------------------------------------

domAbNullC$Study_ID <- as.character(domAbNullC$Study_ID)
domAbNullC$Study_ID <- as.character(domAbNullC$Study_ID)

DomAbNullAndActual<- merge(domAbNullC, DominancModel_Lines, by = "Study_ID")#incuding the slopes values from the actual mixed model fo data
names(DomAbNullAndActual)[1006] <- "ActualSlopes"

DomAbNullAndActual[c(1:5),c(1000:1004)]

#how many studies show a noticable change?
#---------------------------------------------


#selecting studies where the dominant abundance changes more than woudl be expected from the null model
DomAbNullAndActual$DominanceChange <- DomAbNullAndActual$Study_ID #making a spare column to input data 


DomAbNullAndActual$DominanceChange[DomAbNullAndActual$ActualSlopes < DomAbNullAndActual$Q2.5] <- "decrease"
DomAbNullAndActual$DominanceChange[DomAbNullAndActual$ActualSlopes > DomAbNullAndActual$Q97.5] <- "increase"
DomAbNullAndActual$DominanceChange[DomAbNullAndActual$ActualSlopes < DomAbNullAndActual$Q97.5 & DomAbNullAndActual$ActualSlopes > DomAbNullAndActual$Q2.5] <- "noChange" 

barplot(table(DomAbNullAndActual$DominanceChange))

#get actual percentiles for plotting
#----------------------------------------

percentilleData <- data.frame(matrix(NA, length(domAbNull$Study_ID), 2))
names(percentilleData) <- c("Study_ID", "relativeDom")
percentilleData[,1] <- domAbNull$Study_ID

for (p in c(1:length(domAbNull$Study_ID))){

	pCurve <- ecdf(domAbNull[p,2:1001])(DomAbNullAndActual$ActualSlopes[p])
	percentilleData[p, 2] <- pCurve
}

names(DomAbNullAndActual)
merge(percentilleData, DomAbNullAndActual[,c(1, 1005)], by = "Study_ID")

write.csv(percentilleData, "relativeDomChange_Percentiles.csv")

#calculate Z scores : (slope - meanslope)/standard deviation slope 
###-------------------------

DomAbNullAndActual$nullmean <- rowMeans(DomAbNullAndActual[,2:1001])
names(DomAbNullAndActual)
DomAbNullAndActual$SD <- apply(DomAbNullAndActual[,2:1001],1,sd)

DomAbNullAndActual$AssSizeZscores <- (DomAbNullAndActual$ActualSlopes - DomAbNullAndActual$nullmean)/DomAbNullAndActual$SD

#write.csv(DomAbNullAndActual[,c(1, 1008)], "relDom_Zscores.csv")
###making a sloep plot for presentations
#----------------------------------------------
names(metaData)
aDomLines <- coef(DominancModel)$Study_ID
aDomLines$Study <- rownames(aDomLines) #extract the ranom effect levels
lines2 <- merge(aDomLines, metaData[,c(1, 3, 15, 14)], by.x = "Study", by.y = "STUDY_ID", all = FALSE) #add this infor to info on line slopes 
lines <- merge(lines2, domAb[,c(5,7,8,9)], by.x = "Study", by.y = "Study_ID")

head(lines)
lines <-lines[!duplicated(lines), ]
lines$MeanYear <- lines$START_YEAR + lines$LengthYear/2
names(lines)[2] <- "intercept"
names(lines)[3] <- "slope"
head(lines)
nrow(lines)
lines$Study

#work out start and end point of my studies (random effects)  
lines$meanCentedStartYear <- lines$START_YEAR - lines$MeanYear
lines$meanCentedEndYear <- lines$END_YEAR - lines$MeanYear
lines$StartYActuall <- lines$slope * lines$meanCentedStartYear + lines$intercept
lines$EndYActuall <- lines$slope * lines$meanCentedEndYear + lines$intercept
domAb$ActualYear <- domAb$MeanYear + domAb$MeanCentredYear

#work out start and end point of my studies of actuall dates (random effects)  
lines$StartY <- lines$slope * lines$meanCentedStartYear + lines$intercept
lines$EndY <- lines$slope * lines$meanCentedEndYear + lines$intercept

summary(DominancModel)
#manualy calculate the mean slope
endypoint <- -0.044  * 127 + 13.4987
head(domAb)

#plot the same but with points see through and no colour differences 
m3.plot2 <- ggplot(data = domAb, aes(x = ActualYear, y = logDomAb)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plotA2 <- m3.plot2 + theme_classic()+
	theme(text = element_text(size=30)) +
	geom_segment(data =  lines, aes(x = START_YEAR, y = lines$StartY, xend = lines$END_YEAR, yend = lines$EndY),alpha = 2/3, colour = "#333333")+
	geom_segment(aes(x = 1899, y = 13.4987, xend = 2017, yend = endypoint), colour = "red",lwd=1.25) + 
	labs(x="Year",y="Relative abundance (%)")+
	scale_x_continuous(limits = c(1899, 2017))+ 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank()) 


lines[lines$Study == 215,]
sort(as.numeric(lines$Study))

ggplot(data = lines, aes(x = slope)) + 
	geom_histogram()+
	theme_classic()+
	labs(x = "Change in relative dominance", y ="Number of assemblages")+
	theme(text = element_text(size=30)) + 
	geom_vline(xintercept = 0, linetype="dashed", 
                color = "grey", size=1.5)


lines[lines$slope == max(lines$slope),]

#plot mean % dominance in each assembleg for the supplimentary material
#-----------------------------------------------------------

names(domAb)
domAb$percentAb[domAb$Study_ID == "360"]
plotRelAb <- data.frame(domAb %>%
	group_by(Study_ID)%>%
	summarise(meanRelAb = mean(percentAb)))


relHist <- ggplot(data = plotRelAb, aes(meanRelAb))
relHist + geom_histogram()+
	theme_classic()+
	theme(text = element_text(size=30))+
	labs(x = "Mean Relative Dominance (%)", y ="Number of assemblages")+ 
	xlim(5, 100)

length(plotRelAb$meanRelAb[plotRelAb$meanRelAb <21])










