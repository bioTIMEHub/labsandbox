setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
library(lme4)
library(ggplot2)

#
#IS evenness sytematically changing? Evenness calcuated from relative BIOMASS
#------------------------------------------------------------------------------


AbundanceData <- read.csv("evennessAbundanceSimp.csv")
names(AbundanceData)[4] <- "StudyID"
metaData <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data\\metaDataNov.csv")
names(metaData)

head(AbundanceData)

#adding necissary info from metatdata to rarefied data
metaDataRealm <- metaData[,c(1,8, 13, 14, 2, 3)]

AbunanceRarefied <- merge(AbundanceData, metaDataRealm, by.x ="StudyID", by.y = "STUDY_ID" )
head(AbunanceRarefied)
AbunanceRarefied$LengthYears <- AbunanceRarefied$END_YEAR - AbunanceRarefied$START_YEAR

dataRarefied <- AbunanceRarefied

dataRarefied$START_YEAR[dataRarefied$StudyID == 421]
min(dataRarefied$Year[dataRarefied$StudyID == 421])

dataRarefied$START_YEAR[dataRarefied$StudyID == 421] <- 1999

#removing mistakes in taxa list
dataRarefied$TAXA <- as.character(dataRarefied$TAXA)
dataRarefied$TAXA[dataRarefied$TAXA == "Marine invertebrates"] <- "Invertebrates" 
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial plants"] <- "Plants/Algae" 
dataRarefied$TAXA[dataRarefied$TAXA == "Marine plants"] <- "Plants/Algae"  
dataRarefied$TAXA[dataRarefied$TAXA == "Marine invertebrates/plants"] <- "Multipe" 
dataRarefied$TAXA[dataRarefied$TAXA == "Marine Invertebrates"] <- "Invertebrates" 
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial Plants"] <- "Plants/Algae"
dataRarefied$TAXA[dataRarefied$TAXA == "Freshwater invertebrates"] <- "Invertebrates"
dataRarefied$TAXA[dataRarefied$TAXA == "Terrestrial invertebrates"] <- "Invertebrates"
dataRarefied$TAXA <- as.factor(dataRarefied$TAXA)
levels(dataRarefied$TAXA)

#centering each dataset around the mean so that my intercept value will be my mean value as well in the model
dataRarefied$MeanCentredYear <- dataRarefied$Year

for (i in unique(dataRarefied$Study)){
	studyi <- dataRarefied[dataRarefied$Study == i,]	# select data from a singel study
	meani <- mean(studyi$Year)	#calculate mean year
	dataRarefied$MeanYear[dataRarefied$Study== i] <- meani
	meanCentredYeari <- studyi$Year - meani # subtracting the mean from each year 
	dataRarefied$MeanCentredYear [dataRarefied$Study == i] <- meanCentredYeari		#putting data into the dataset
}


dataRarefied$Study <- as.factor(dataRarefied$Study) # ensure study is a factor not a numeric
dataRarefied$meanCentredStartYear <- dataRarefied$START_YEAR - dataRarefied$MeanYear 
dataRarefied$meanCentredEndYear <- dataRarefied$END_YEAR - dataRarefied$MeanYear 

dataRarefied$StudyID[dataRarefied$meanCentredStartYear == min(dataRarefied$meanCentredStartYear)]


dataRarefied$Study <- as.factor(dataRarefied$Study) # ensure study is a factor not a numeric



####log transformed N against mean centred year with varying slope due to study and length of study #####
#-------------------------------------------------------

m3 <- lmer(Evenness~ MeanCentredYear + (1 + MeanCentredYear|Study), data=dataRarefied)
summary(m3)

#check for normal distribution
sresid <- resid(m3, trpe = "pearson")
hist(sresid)	#looks normal

#residuals vs fitted values
fits <- fitted(m3)
plot(sresid ~ fits)	#no problem

#residuals vs each x variable
plot(sresid ~ dataRarefied$MeanCentredYear) # no pattern that I can see
names(dataRarefied)

#plot this relationship
#------------------------------------

names(dataRarefied)

m3_Lines <- coef(m3)$Study[c(1:2)] # get slope and intercept of each line in model 3
names(m3_Lines) <- c("intercept", "slope")

m3_Lines$Study <- rownames(m3_Lines) #extract the ranom effect levels
lines <- merge(m3_Lines, dataRarefied[,c(1, 5, 6, 14, 15)], by.x = "Study", by.y = "StudyID", all = FALSE) #add this infor to info on line slopes 
head(lines)
lines <- lines[!duplicated(lines), 

lines$START_YEAR[lines$Study == 421]

names(dataRarefied)
min(dataRarefied$MeanCentredYear)

#work out start and end point of my studies (random effects)  
lines$StartY <- lines$slope * lines$meanCentredStartYear + lines$intercept
lines$EndY <- lines$slope * lines$meanCentredEndYear + lines$intercept

m3.plot <- ggplot(data = dataRarefied, aes(x = MeanCentredYear, y = Evenness)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plot + geom_point() + 
	geom_segment(data =  lines, aes(x = meanCentredStartYear, y = StartY, xend = meanCentredEndYear, yend = EndY))


#plot the same but with points see through and no colour differences 
m3.plot <- ggplot(data = dataRarefied, aes(x = MeanCentredYear, y = AbudnaceLog)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plotA <- m3.plot + theme_classic()+
	theme(text = element_text(size=30)) +
	geom_segment(data =  lines, aes(x = meanCentredStartYear, y = StartY, xend = meanCentredEndYear, yend = EndY),alpha = 2/3)+
	labs(x="Mean Centred Year",y="Evenness")  

ggsave("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\EvennesChangemcy.png", bg = "transparent")

###Making a histogram of slopes of change
#---------------------------------------------

hist(lines$slope)

histPlot <- ggplot(aes(x = slope), data = lines)
histPlot <- histPlot + geom_histogram() +
	theme_classic() +
	labs(x = "Rate of Change of Evenness", y = "Number of Assemblages")+
	theme(text = element_text(size=25)) + 
	theme(
  	  panel.grid.minor = element_blank(), 
  	  panel.grid.major = element_blank(),
  	  panel.background = element_blank(),
  	  plot.background = element_blank())+
	geom_vline(aes(xintercept= 0),
            color="grey", linetype="dashed", size=1)
	#+
	#geom_vline(aes(xintercept=0.002730),
      #      color="red", linetype="dashed", size=1)

write.csv(lines, "evennessChangeSimpAb.csv")


#Making a nice line plot
#-----------------------------

MyGreen <- rgb(72/225,203/225,61/225)

names(dataRarefied)
names(m3_Lines)

m3_Lines$Study <- rownames(m3_Lines) #extract the ranom effect levels
lines <- merge(m3_Lines, dataRarefied[,c(1, 5, 6, 7, 14, 15)], by.x = "Study", by.y = "StudyID", all = FALSE) #add this infor to info on line slopes 
lines <- lines[!duplicated(lines), ]

#work out start and end point of my studies (random effects)  
lines$StartYActuall <- lines$slope * lines$START_YEAR + lines$intercept
lines$EndYActuall <- lines$slope * lines$END_YEAR + lines$intercept
dataRarefied$ActualYear <- dataRarefied$MeanYear + dataRarefied$MeanCentredYear

#work out start and end point of my studies of actuall dates (random effects)  
lines$StartY <- lines$slope * lines$meanCentredStartYear + lines$intercept
lines$EndY <- lines$slope * lines$meanCentredEndYear + lines$intercept


#manualy calculate the mean slope
endypoint <- -0.0004410 * 127 + 0.6044398


m3.plot2 <- ggplot(data = dataRarefied, aes(x = ActualYear, y = Percent)) # plot of model with varying intercept but not slope, but only main slope shown
m3.plotA2 <- m3.plot2 + theme_classic()+
	theme(text = element_text(size=30)) + 
	labs(x="Year",y="Evenness")+
	scale_x_continuous(limits = c(1899, 2017))+
	geom_segment(data =  lines, aes(x =  lines$START_YEAR, y = lines$StartY, xend =  lines$END_YEAR, yend = lines$EndY),alpha = 2/3, colour = "#333333")+
	geom_segment(aes(x = 1899, y = 0.6044398, xend = 2017, yend = endypoint), colour = MyGreen,lwd=1.25)#









