#script for creating a plot for my dominance paper. I repeat the analysis 
#undertaken in the two seperate relative and abosolute doimunance scripts
#that dont involve teh null model, and then put plots 
#together in a jpeg

setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\NullModel_correctsubset")
AbundanceDataAll <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\data preperation\\AbudnaceOnceRarefy.csv")
names(AbundanceDataAll)

metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Data_Oct2017\\bioTIMEmetadataJune.csv")


library(matrixStats)
library(lme4)
library(vegan)
library(paleotree)
library(tidyr)
library(binhf)
library(reshape2)
library(dplyr)
library(ggplot2)
library(gridExtra)

#running mixed model on actual dominance data 
#--------------------------------------------

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

#calculate mean centred year for the model
#AbundanceDataAll2$MeanCentredYear <- AbundanceDataAll2$Year - AbundanceDataAll2$MeanYear

#head(AbundanceDataAll)

#count length of studies
AbSelect <- data.frame(AbundanceDataAll2 %>%
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
rarefiednopaS <- rarefiednopa[!rarefiednopa$Study_ID %in% removeSTudies,]


head(AbSelect)

#get the dominant species each year and output its abundance and percentage abundance
domAb <- data.frame(rarefiednopaS %>%
		group_by(Study_ID, MeanCentredYear) %>%
		mutate(abSum = sum(Abudance)) %>% #get overall abunacne of all individulas that year/study
		filter(Abudance == max(Abudance)[1]) %>% #select the largest value (if there is more than one tied then just chose the first value)
		filter(1:n() == 1) %>% #make sure  there is only oen value!
		mutate(logDomAb = log2(Abudance + 1), #calculate log abundance
			percentAb = (Abudance/abSum * 100)) #calculate percentage abundance
	)

head(domAb)


#run model for relative dominance
#---------------------------------- 
relDominancModel <- lmer(percentAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAb)
summary(relDominancModel )

relDominancModel_Lines <- coef(relDominancModel)$Study
relDominancModel_Lines$Study_ID <- rownames(relDominancModel_Lines)
names(relDominancModel_Lines) <- c("rel_intercept", "rel_Slope", "Study_ID")

#run model for absolute dominance
#---------------------------------- 
AbDominancModel <- lmer(logDomAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAb)
summary(AbDominancModel )

AbDominancModel_Lines <- coef(AbDominancModel)$Study
AbDominancModel_Lines$Study_ID <- rownames(AbDominancModel_Lines)
names(AbDominancModel_Lines) <- c("ab_intercept", "ab_Slope", "Study_ID")


#merge both values
linesAll <- merge(relDominancModel_Lines, AbDominancModel_Lines, by = "Study_ID")

###making a sloep plot for presentations
#----------------------------------------------

names(metaData)

lines2 <- merge(linesAll, metaData[,c(1, 3, 15, 14)], by.x = "Study_ID", by.y = "STUDY_ID", all = FALSE) #add this infor to info on line slopes 
lines <- merge(lines2, domAb[,c(5,7,8,9)], by.x = "Study_ID")

head(lines)
lines <-lines[!duplicated(lines), ]
lines$MeanYear <- lines$StartYear + lines$nyear/2
head(lines)
nrow(lines)
lines$Study

#Plot Absolute Domiance 
#------------------------

#work out start and end point of my studies (random effects)  
lines$meanCentedStartYear <- lines$StartYear - lines$MeanYear
lines$meanCentedEndYear <- lines$endYear - lines$MeanYear
domAb$ActualYear <- domAb$MeanYear + domAb$MeanCentredYear

#work out start and end point of my studies of actuall dates (random effects)  
lines$ab_StartY <- lines$ab_Slope * lines$meanCentedStartYear + lines$ab_intercept
lines$ab_EndY <- lines$ab_Slope * lines$meanCentedEndYear + lines$ab_intercept

summary(DominancModel)

ab_endypoint <- 0.005816  * 127 + 9.117

#plot the same but with points see through and no colour differences 
AbDomPlot <- ggplot(data = domAb, aes(x = ActualYear, y = logDomAb)) # plot of model with varying intercept but not slope, but only main slope shown
AbDomPlot <- AbDomPlot + theme_classic()+
	theme(text = element_text(size=20)) +
	geom_segment(data = lines, aes(x = StartYear, y = lines$ab_StartY, xend = lines$endYear, yend = lines$ab_EndY),alpha = 2/3, colour = "#333333")+
	geom_segment(aes(x = 1899, y = 9.117, xend = 2017, yend = ab_endypoint), colour = "blue",lwd=1.25) + 
	labs(x="",y="Log Abundance")+
	scale_x_continuous(limits = c(1899, 2017))+ 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank())+
  	geom_text(aes(x = 1900, y = 26), label = "A)", size = 6)
 


#plot relative abundance
#------------------------

#work out start and end point of my studies of actuall dates (random effects)  
lines$rel_StartY <- lines$rel_Slope * lines$meanCentedStartYear + lines$rel_intercept
lines$rel_EndY <- lines$rel_Slope * lines$meanCentedEndYear + lines$rel_intercept

rel_endypoint <- -0.044  * 127 + 13.4987

#plot the same but with points see through and no colour differences 
RelDomPlot <- ggplot(data = domAb, aes(x = ActualYear, y = percentAb)) # plot of model with varying intercept but not slope, but only main slope shown
RelDomPlot <- RelDomPlot + theme_classic()+
	theme(text = element_text(size=20)) +
	geom_segment(data =  lines, aes(x = StartYear, y = lines$rel_StartY, xend = lines$endYear, yend = lines$rel_EndY),alpha = 2/3, colour = "#333333")+
	geom_segment(aes(x = 1899, y = 13.4987, xend = 2017, yend = rel_endypoint), colour = "red",lwd=1.25) + 
	labs(x="Year",y="Relative abundance (%)")+
	scale_x_continuous(limits = c(1899, 2017))+ 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank()) +
  	geom_text(aes(x = 1900, y = 76), label = "B)", size = 6)


bothPlot <- grid.arrange(AbDomPlot, RelDomPlot, ncol = 1)

ggsave("LinePlot.jpeg", plot = bothPlot, device = "jpeg", width = 13, height = 20, units = c("cm"))









