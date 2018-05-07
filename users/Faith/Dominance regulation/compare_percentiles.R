setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\NullModel_correctsubset")

#code to compare the perentiles of each of dominance change and assemblage chaneg measurements together.
#data was writen in the "AssemblageSize.r, ActualDominance.r and PercentDominance.r sripts. Each number
#is what percentile of the null distribution the actual sloep of change is. 

library(ggplot2)
library(dplyr)

AbundanceDataAll <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\data preperation\\AbudnaceOnceRarefy.csv")
names(AbundanceDataAll)

#importing perecntiles data and setting up data frame
#------------------------------------------------------
relPe <- read.csv("relativeDomChange_Percentiles.csv")
AssPe <- read.csv("AssSizeChange_Percentiles.csv")
AcrDPe <- read.csv("ActualDomChange_Percentiles.csv")

relPe$X <- NULL
AssPe$X <- NULL
AcrDPe$X <- NULL

relaAss <- merge(relPe, AssPe, by = "Study_ID")
AllPercentiles <- merge(relaAss, AcrDPe, by = "Study_ID")
names(AllPercentiles)

#importing z scores and setting up data frame 
#--------------------------------------------
SRZscore <- read.csv("SpecRich_Zscores.csv")
AssSizeZscore <- read.csv("AssSize_Zscores.csv")
RelDomZscore <- read.csv("RelDom_Zscores.csv")
ActDomZscore <- read.csv("ActDom_Zscores.csv")

SRZscore$X <- NULL
AssSizeZscore$X <- NULL
RelDomZscore$X <- NULL
ActDomZscore$X <- NULL
names(RelDomZscore)[2] <- "relDomZscores"

#merge info together
SRass<- merge(SRZscore, AssSizeZscore, by = "Study_ID")
SRassRel<- merge(RelDomZscore, SRass, by = "Study_ID")
AllZscores<- merge(SRassRel, ActDomZscore, by = "Study_ID")


#plotting percentile data against each other
#--------------------------------------------

plot(AllPercentiles$relativeDom ~ AllPercentiles$ActualDom, 
	xlab = "Actual Dominant Abundance Change", ylab = "Relative Dominant Abundance Change")
cor(AllPercentiles[,c(2,4)], method = "spearman")

plot(AllPercentiles$AssSize ~ AllPercentiles$relativeDom, 
	xlab = "Relative Dominant Abundance", ylab = "Assemblage size change")
cor(AllPercentiles[,c(2,3)], method = "spearman")


plot(AllPercentiles$AssSize ~ AllPercentiles$ActualDom, 
	xlab = "Actual Dominant Abundance Change", ylab = "Assemblage size change")
cor(AllPercentiles[,c(3,4)], method = "spearman")


###add extra info on species richnes and assembalges 
#-----------------------------------------------

#count length of studies
AbSelect <- data.frame(AbundanceDataAll %>%
	group_by(Study_ID)%>%
	mutate(LengthYear = n_distinct(Year))%>%
	filter(LengthYear > 9)
)

head(AbSelect)

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



#get the dominant species each year and output its abundance and percentage abundance
domAb <- data.frame(rarefiednopaS %>%
		group_by(Study_ID, Year) %>%
		mutate(abSum = sum(Abudance),
			SpeciesRichness = n_distinct(Species_Identity)) %>% #get overall abunacne of all individulas that year/study
		filter(Abudance == max(Abudance)[1]) %>% #select the largest value (if there is more than one tied then just chose the first value)
		filter(1:n() == 1) %>% #make sure  there is only oen value!
		mutate(logDomAb = log2(Abudance + 1), #calculate log abundance
			percentAb = (Abudance/abSum * 100),
			logAssSize = log2(abSum + 1)) #calculate percentage abundance
	)

head(domAb)
plotData <- data.frame(domAb %>%
	group_by(Study_ID) %>%
	summarise(meanSR = mean(SpeciesRichness), meanAssSize = mean(abSum),
		nYear = mean(LengthYear)))

#for percentiles
#-----------------
#combine data for percentiles
AllPlotData <- merge(AllPercentiles, plotData, by = "Study_ID")

#add columns saying whether the value was above or below the 95% percentiles

AllPlotData$relativeChange <- "NoChange"
AllPlotData$AbDomChange <- "NoChange"

AllPlotData$relativeChange[which(AllPlotData$relativeDom > 0.975)] <- "increase"
AllPlotData$relativeChange[which(AllPlotData$relativeDom < 0.025)] <- "decrease"
AllPlotData$AbDomChange[which(AllPlotData$ActualDom > 0.975)] <- "increase"
AllPlotData$AbDomChange[which(AllPlotData$ActualDom < 0.025)] <- "decrease"

#create plots for percentiles
#---------------------------------------------

head(AllPlotData)
AbDom <- ggplot(data = AllPlotData, aes(x = ActualDom, y = nYear, size = meanSR, colour = AbDomChange))
AbDom + geom_point(alpha=.4)+
	theme_bw() + 
	scale_color_manual(values=c("Purple", "Green", "black"))+
	labs(x = "Absolute Dominance Change Percentile", y = "Number of Sampling Years")+
	theme(text = element_text(size=20))

head(AllPlotData)
relDom <- ggplot(data = AllPlotData, aes(x = relativeDom, y = nYear, size = meanSR, colour = relativeChange))
relDom + geom_point(alpha=.4)+
	theme_bw() + 
	scale_color_manual(values=c("Purple", "Green", "black"))+
	labs(x = "Relative Dominance Change Percentile", y = "Number of Sampling Years")+
	theme(text = element_text(size=20))

#for Z scores
#-----------------
#combine data for percentiles and z scores so i can see what studies had significan changes 
AllZscoresPlot <- merge(AllZscores, AllPlotData, by = "Study_ID")

#create plots for Z scores 
#---------------------------------------------

head(AllZscores)
AbDomz <- ggplot(data = AllZscoresPlot, aes(x = ActDomZscores, y = nYear, size = meanSR, colour = AbDomChange))
AbDomz + geom_point(alpha=.4)+
	theme_bw() + 
	scale_color_manual(values=c("Purple", "Green", "black"))+
	labs(x = "Absolute Dominance Change Z score", y = "Number of Sampling Years")+
	theme(text = element_text(size=20))

head(AllPlotData)
relDomz <- ggplot(data = AllZscoresPlot, aes(x = relDomZscores, y = nYear, size = meanSR, colour = relativeChange))
relDomz + geom_point(alpha=.4)+
	theme_bw() + 
	scale_color_manual(values=c("Purple", "Green", "black"))+
	labs(x = "Relative Dominance Z score", y = "Number of Sampling Years")+
	theme(text = element_text(size=20))





















