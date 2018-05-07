#thsi code takes teh correct subset of data (not including presence/absence data) and sees if there are any interesting 
#patterns in splitting by taxa, realm or experiment 


setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Null Model")
domAbNull <- read.csv("DominantAbundanceNull.csv")
domAbNull[,2]
head(domAbNull)

AbundanceDataAll <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\data preperation\\AbudnaceOnceRarefy.csv")
names(AbundanceDataAll)
metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Data_Oct2017\\bioTIMEmetadataJune.csv")
names(metaData)
table(metaData$TAXA)
#compare the mixed model of actual abundance of dominant species with the null model results
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
library(lmerTest)

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
rarefiednopaS <- rarefiednopa[!rarefiednopa$Study_ID %in% removeSTudies,]


#getting correct taxa labels 
#---------------------------------
metaData$TAXA2 <- as.character(metaData$TAXA)
unique(metaData$TAXA2)
metaData$TAXA2[metaData$TAXA2 == "Marine invertebrates"] <- "Invertebrates" 
metaData$TAXA2[metaData$TAXA2 == "Terrestrial plants"] <- "Plants/Algae" 
metaData$TAXA2[metaData$TAXA2 == "Marine plants"] <- "Plants/Algae"  
metaData$TAXA2[metaData$TAXA2 == "Marine invertebrates/plants"] <- "Multipe" 
metaData$TAXA2[metaData$TAXA2 == "Marine Invertebrates"] <- "Invertebrates" 
metaData$TAXA2[metaData$TAXA2 == "Terrestrial Plants"] <- "Plants/Algae"
metaData$TAXA2[metaData$TAXA2 == "Freshwater invertebrates"] <- "Invertebrates" 
metaData$TAXA2[metaData$TAXA2 == "Terrestrial invertebrates"] <- "Invertebrates" 
metaData$TAXA2[metaData$TAXA2 == "All"] <- "Multipe" 
metaData$TAXA2[metaData$TAXA2 == "Freshwater plants"] <- "Plants/Algae"

#running mixed model on actual dominance data 
#--------------------------------------------

table(metaData$TAXA2)

head(rarefiednopa)
rarefiednopa$Study_ID[rarefiednopa$Abudance == min(rarefiednopa$Abudance)]

#get the dominant species each year and output its abundance and percentage abundance
domAb <- data.frame(rarefiednopaS %>%
		group_by(Study_ID, MeanCentredYear) %>%
		mutate(abSum = sum(Abudance),
			SpeciesRichness = n_distinct(Species_Identity)) %>% #get overall abunacne of all individulas that year/study
		filter(Abudance == max(Abudance)[1]) %>% #select the largest value (if there is more than one tied then just chose the first value)
		filter(1:n() == 1) %>% #make sure  there is only oen value!
		mutate(logDomAb = log2(Abudance + 1), #calculate log abundance
			percentAb = (Abudance/abSum * 100),
			logAssSize = log2(abSum + 1)) #calculate percentage abundance
	)

head(domAb)
domAb <- domAb[!domAb$logDomAb == 0,]#remove study 360 wich has biomass data not abundance data


#add realm and taxa data
domAbMeta <- merge(domAb, metaData[,c(1, 3, 28)], by.x = "Study_ID", by.y = "STUDY_ID")
names( metaData)

#run model of SpeciesRichness change
SpeciesRichnessModel <- lmer(SpeciesRichness ~ MeanCentredYear  + (1 + MeanCentredYear|Study_ID), data=domAb)
summary(SpeciesRichnessModel)
SpeciesRichnessModel_Lines <- coef(SpeciesRichnessModel)$Study[2]
SpeciesRichnessModel_Lines$Study_ID <- rownames(SpeciesRichnessModel_Lines)
names(SpeciesRichnessModel_Lines)[1] <- "SpeciesRichnessChange"

#run model of assemblage size change
AssSizeModel <- lmer(logAssSize ~ MeanCentredYear  + (1 + MeanCentredYear|Study_ID), data=domAb)
summary(AssSizeModel)
AssSizeModel_Lines <- coef(AssSizeModel)$Study[2]
AssSizeModel_Lines$Study_ID <- rownames(AssSizeModel_Lines)
names(AssSizeModel_Lines)[1] <- "AssSizeChange"


#add SR and ass size slopes into main dataset
bothLines <- merge(AssSizeModel_Lines, SpeciesRichnessModel_Lines, by = "Study_ID")
domAbAss <- merge(domAbMeta, bothLines, by = "Study_ID")
tail(domAbAss)

#ABSOULTE DOMIANCE

#absolute dominance with assembalge size chaneg and REALM
DominancModelR <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + REALM + REALM:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModelR)

DominancModelR2 <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + REALM:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModelR2)
anova(DominancModelR,DominancModel2)

#absolute dominance with assembalge size chaneg and TAXA
DominancModelT <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + TAXA2 + TAXA2:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModelT)

# absolute dominance against year WITHOUT an interaction effect of assemblage size change or specioes richness
DominancModel2 <- lmer(logDomAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModel2)

# absolute dominance against year with an interaction effect of assemblage size change but NOT species richness
DominancModel3 <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModel3)

DomSlope <- coef(DominancModel3)$Study_ID
DomSlope$Study_ID <- rownames(DomSlope)

# absolute dominance against year with an interaction effect of assemblage size change and species richness
DominancModel <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + SpeciesRichnessChange + 
				SpeciesRichnessChange:MeanCentredYear+(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModel)
anova(DominancModel)

# absolute dominance against year with an interaction effect of assemblage size change and species richness, and these interacting together
DominancModel4 <- lmer(logDomAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + SpeciesRichnessChange + 
				SpeciesRichnessChange:MeanCentredYear + SpeciesRichnessChange:AssSizeChange +
				(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(DominancModel4)

#compare model fits
anova(DominancModel2, DominancModel3)#the model with assemblge size is definately better than the one without 
anova(DominancModel2, DominancModel)
anova(DominancModel4, DominancModel)

#RELATIVE DOMAINANCE

#dominance against year with NO interaction effect of assemblage size change or species richness
relDomModel3 <- lmer(percentAb ~ MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModel3)
anova(relDomModel3)
RelDomSlope <- coef(relDomModel3)$Study_ID
RelDomSlope$Study_ID <- rownames(RelDomSlope)

#relative dominance with assembalge size chaneg and REALM
relDomModelR <- lmer(percentAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + REALM + REALM:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModelR)

relDomModelR2 <- lmer(percentAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + REALM:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModelR2)
anova(relDomModelR,relDomModelR2)

relDomModelT <- lmer(percentAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + TAXA2 + TAXA2:MeanCentredYear + (1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModelT)


#dominance against year with NO interaction effect of assemblage size change or species richness, but assSize as a fixed effect
relDomModel5 <- lmer(percentAb ~ MeanCentredYear + AssSizeChange +(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModel5)
anova(relDomModel5)

#dominance against year with an interaction effect of assemblage size change 
relDomModel2 <- lmer(percentAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModel2)
anova(relDomModel2)

#relative dominance against year with an interaction effect of assemblage size change and species richness
relDomModel <- lmer(percentAb ~ MeanCentredYear + AssSizeChange + AssSizeChange:MeanCentredYear + SpeciesRichnessChange + 
				SpeciesRichnessChange:MeanCentredYear+(1 + MeanCentredYear|Study_ID), data=domAbAss)
summary(relDomModel)
anova(relDomModel)

#compare model fits
anova(relDomModel3, relDomModel2)
anova(relDomModel3, relDomModel)
anova(relDomModel3, relDomModel5)


#plotting boxplots of slopes slpit by REALM and TAXA
#-----------------------------------------------

Ab_lines <- coef(DominancModel2)$Study
names(Ab_lines) <- c("InterceptAb", "Slope_ab")
Ab_lines$Study_ID <- rownames(Ab_lines) 

Rel_lines <- coef(relDomModel5)$Study
names(Rel_lines) <- c("InterceptRel", "Slope_Rel")
Rel_lines$Study_ID <- rownames(Rel_lines)


Lines <- merge(Ab_lines, Rel_lines, by = "Study_ID")
PlotLines <- merge(Lines, metaData[,c(1, 3, 28)], by.x = "Study_ID", by.y = "STUDY_ID")

#plotting boxplots fo realm
abRealmPlot <- ggplot(aes(x = REALM, y = Slope_ab), data = PlotLines)
abRealmPlot + geom_boxplot()+
	geom_hline(yintercept = 0, linetype = 2, colour = "grey")+ 
	geom_boxplot() +
	theme_classic() +
	labs(x = "Realm", y ="Rate of Change of Absolute Dominance")+
	theme(text = element_text(size=25)) 

relRealmPlot <- ggplot(aes(x = REALM, y = Slope_Rel), data = PlotLines)
relRealmPlot + geom_boxplot()+
	geom_hline(yintercept = 0, linetype = 2, colour = "grey") + 
	geom_boxplot() +
	theme_classic() +
	labs(x = "Realm", y ="Rate of Change of Relative Dominance")+
	theme(text = element_text(size=25)) 


#plotting boxplots fo taxa
abTaxaPlot <- ggplot(aes(x = TAXA2, y = Slope_ab), data = PlotLines)
abTaxaPlot + geom_boxplot()+
	geom_hline(yintercept = 0, linetype = 2, colour = "grey") + 
	geom_boxplot() + 
	theme_classic() +
	labs(x = "", y ="Rate of Change of Absolute Dominance")+
	theme(text = element_text(size=25)) 

relTaxaPlot <- ggplot(aes(x = TAXA2, y = Slope_Rel), data = PlotLines)
relTaxaPlot + geom_boxplot()+ 
	geom_hline(yintercept = 0, linetype = 2, colour = "grey") + 
	geom_boxplot() +
	theme_classic() +
	labs(x = "Taxa", y ="Rate of Change of Relative Dominance")+
	theme(text = element_text(size=25)) 


####Running models without mixed taxa, to check this doesnt chaneg the answer
#-----------------------------------------------------------------------------------

names(domAbAss)
noMultiple <- domAbAss[!domAbAss$TAXA2 == "Multipe", ]

# absolute dominance against year 
noMultipleModel <- lmer(logDomAb ~ MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=noMultiple)
summary(noMultipleModel)

# relative dominance against year 
noMultipleModelpercent <- lmer(percentAb ~ MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=noMultiple)
summary(noMultipleModelpercent)

####Running models without mixed taxa or benthic communities, to check this doesnt chaneg the answer
#-----------------------------------------------------------------------------------
names(domAbAss)
noMultipleBenth <- domAbAss[!domAbAss$TAXA2 == "Multipe"|!domAbAss$TAXA2 == "Benthos", ]

unique(noMultipleBenth$TAXA2)

# absolute dominance against year 
noMultipleBenthModel <- lmer(logDomAb ~ MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=noMultipleBenth)
summary(noMultipleBenthModel)

# relative dominance against year 
noMultipleBenthModelpercent <- lmer(percentAb ~ MeanCentredYear +(1 + MeanCentredYear|Study_ID), data=noMultipleBenth)
summary(noMultipleBenthModelpercent)

####Running models without the experimental datasets
#-----------------------------------------------------------------------------------

experimental <- c("44", "59", "214", "221", "248", "300", "313", "336") 
names(domAbAss)
noExperiment <- domAbAss[!domAbAss$Study_ID %in% experimental, ]

# absolute dominance against year 
noExperimentModel <- lmer(logDomAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data = noExperiment)
summary(noExperimentModel)

# absolute dominance against year 
noExperimentModelpercent <- lmer(percentAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data = noExperiment)
summary(noExperimentModelpercent)

####Running models without the experimental study 248
#-----------------------------------------------------------------------------------

noExperiment248 <- domAbAss[!domAbAss$Study_ID == "248", ]

# absolute dominance against year 
noExperiment248Model <- lmer(logDomAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data = noExperiment248)
summary(noExperiment248Model)

# absolute dominance against year 
noExperiment248Modelpercent <- lmer(percentAb ~ MeanCentredYear + (1 + MeanCentredYear|Study_ID), data = noExperiment248)
summary(noExperiment248Modelpercent)













