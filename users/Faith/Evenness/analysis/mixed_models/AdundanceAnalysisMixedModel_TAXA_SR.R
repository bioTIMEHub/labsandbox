setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Mixed model")
library(lme4)
library(lsmeans)
library(ggplot2)
library(MuMIn)

#### ABUNDANCE ANALYSIS With REALM####
#---------------------------------------------

dataRarefied <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data\\rarefiedAbundanceNoSmall.csv")
head(dataRarefied)
names(dataRarefied)
metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data\\metaDataNov.csv")
names(metaData)

#adding necissary info from metatdata to rarefied data
metaDataRealm <- metaData[,c(1,2,3,6, 8,13, 14, 15)]
dataRarefied <- merge(dataRarefied, metaDataRealm, by.x ="Study", by.y = "STUDY_ID" )
head(dataRarefied)
dataRarefied$LengthYears <- dataRarefied$END_YEAR - dataRarefied$START_YEAR
problemStudy <- as.character(dataRarefied$Study[dataRarefied$LengthYears > 1000]) # sorting out a problem study 
#
dataRarefied[dataRarefied$Study == problemStudy,19] <- 1999 

dataRarefied$LengthYears <- dataRarefied$END_YEAR - dataRarefied$START_YEAR


#dataRarefied$REALM.y <- dataRarefied$REALM

#removing mistakes in taxa list and group taxa without realm
unique(dataRarefied$TAXA)
dataRarefied$TAXA <- as.character(dataRarefied$TAXA)
dataRarefied$TAXA[dataRarefied$TAXA == "Marine Invertebrates"] <- "Invertebrates" 
dataRarefied$TAXA[dataRarefied$TAXA == "Marine plants"] <- "Plants/Algae"
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial Invertebrates"] <- "Invertebrates"
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial Plants"] <- "Plants/Algae"
dataRarefied$TAXA[dataRarefied$TAXA == "Marine invertebrates"] <- "Invertebrates"
dataRarefied$TAXA[dataRarefied$TAXA == "Marine invertebrates/plants"] <- "Multipe"
dataRarefied$TAXA[dataRarefied$TAXA == "All"] <- "Multipe"
dataRarefied$TAXA[dataRarefied$TAXA == "Marine Plants"] <- "Plants/Algae"
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial invertebrates"] <- "Invertebrates"
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial plants"] <- "Plants/Algae"
dataRarefied$TAXA[dataRarefied$TAXA == "Freshwater invertebrates"] <- "Invertebrates" 


dataRarefied$TAXA <- as.factor(dataRarefied$TAXA)
levels(dataRarefied$TAXA)


# random effects group = dataRarefied$Study
# y = dataRarefied$N
# x = dataRarefied$year 

dataRarefied$LogN <- log1p(dataRarefied$N)

#centering each dataset around the mean so that my intercept value will be my mean value as well in the model
dataRarefied$MeanCentredYear <- dataRarefied$year
dataRarefied$MeanYear <- dataRarefied$year

for (i in unique(dataRarefied$Study)){
	studyi <- dataRarefied[dataRarefied$Study == i,]	# select data from a singel study
	meani <- mean(studyi$year)	#calculate mean year
	dataRarefied$MeanYear[dataRarefied$Study== i] <- meani
	meanCentredYeari <- studyi$year - meani # subtracting the mean from each year 
	dataRarefied$MeanCentredYear [dataRarefied$Study == i] <- meanCentredYeari		#putting data into the dataset
}

head(dataRarefied)

dataRarefied$Study <- as.factor(dataRarefied$Study) # ensure study is a factor not a numeric
dataRarefied$meanCentredStartYear <- dataRarefied$START_YEAR - dataRarefied$MeanYear 
dataRarefied$meanCentredEndYear <- dataRarefied$END_YEAR - dataRarefied$MeanYear 

dataRarefied$LogSR <- log(dataRarefied$S)

##### calculate speceis richness change
#--------------------------------------------------------
SRm3 <- lmer(LogSR  ~ MeanCentredYear + TAXA + (1 + MeanCentredYear|Study), data=dataRarefied)
summary(SRm3)

sresidSRm3 <- resid(SRm3, trpe = "pearson")
hist(sresidSRm3)	#looks normal


m3_LinesSR <- coef(SRm3)$Study[c(1:2)] # get slope and intercept of each line in model 3
names(m3_LinesSR) <- c("SRintercept", "SRslope")

####log transformed N against mean centred year with varying slope - realm and taxa #####
#-------------------------------------------------------

#model without interaction effect
m3 <- lmer(LogN ~ MeanCentredYear + TAXA + (1 + MeanCentredYear|Study), data=dataRarefied)
summary(m3)


m3_Lines <- coef(m3)$Study[c(1:2)] # get slope and intercept of each line in model 3
names(m3_Lines) <- c("intercept", "slope")

CompareData <- cbind(m3_LinesSR, m3_Lines)
CompareData$Study_ID <- rownames(CompareData)

plot(CompareData$SRslope, CompareData$slope, xlab = "Log Species Richness Change", ylab = "Log Assemblage Size change")
abline(lm(CompareData$slope ~ CompareData$SRslope)) 

cor(CompareData$SRslope, CompareData$slope, method = "spearman")
#write.csv(CompareData, "C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data\\SpeciesRichnessChange.csv")


###Species richness plot


#plot model with varying slopes and intercepts

m3_Lines <- coef(SRm3)$Study[c(1:2)] # get slope and intercept of each line in model 3
names(m3_Lines) <- c("intercept", "slope")



names(dataRarefied)

m3_Lines$Study <- rownames(m3_Lines) #extract the ranom effect levels
lines <- merge(m3_Lines, dataRarefied[,c(1, 20, 21, 27, 28)], by = "Study", all = FALSE) #add this infor to info on line slopes 
lines <- lines[!duplicated(lines), ]

names(dataRarefied)
dataRarefied[dataRarefied$Study == 421, ]
min(dataRarefied$MeanCentredYear)


#work out start and end point of my studies (random effects)  
lines$StartY <- lines$slope * lines$meanCentredStartYear + lines$intercept
lines$EndY <- lines$slope * lines$meanCentredEndYear + lines$intercept

m3.plot <- ggplot(data = dataRarefied, aes(x = MeanCentredYear, y = LogSR)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plot + geom_point() + 
	geom_segment(data =  lines, aes(x = meanCentredStartYear, y = StartY, xend = meanCentredEndYear, yend = EndY))


#plot the same but with points see through and no colour differences 
m3.plot <- ggplot(data = dataRarefied, aes(x = MeanCentredYear, y = LogSR)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plotA <- m3.plot + theme_classic()+
	theme(text = element_text(size=30)) +
	geom_segment(data =  lines, aes(x = meanCentredStartYear, y = StartY, xend = meanCentredEndYear, yend = EndY),alpha = 2/3) + 
	labs(x="Mean Centred Year",y="Log2 Species Richness")


ggsave("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\RichnessChangemcy.png", bg = "transparent")

###Making a histogram of slopes of change
#---------------------------------------------

hist(lines$slope)

histPlot <- ggplot(aes(x = slope), data = lines)
histPlot <- histPlot + geom_histogram() +
	theme_classic() +
	labs(x = "Rate of Change of Species Richness", y = "Number of Studies")+
	theme(text = element_text(size=25)) + 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank())+
	geom_vline(aes(xintercept= 0.0),
            color="grey", linetype="dashed", size=1)







