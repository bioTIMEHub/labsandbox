setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
AbData <- read.csv("AbudnaceOnceRarefy.csv")
library(vegan)
library(reshape2)

#code to calculate evenness change
#----------------------------------------------------------------

head(AbData)
names(AbData)[2] <- "abundance"
#select data to practice on
AbData10 <- AbData[AbData$Study_ID == 10,]
AbData10_y1 <- AbData10[AbData10$Year == AbData10$Year[1],]

#sum abundance for each species 
evennessData <- aggregate(AbData10_y1$abundance, list(AbData10_y1$Species_Identity), sum)

H <- diversity(evennessData[,2], index = "simpson")#Simpson index
J <- H/specnumber(evennessData[,2])#calculate Pileu'e evenness 

#calculate evenness each year in dataset 10 using a loop

EvenData10 <- data.frame(unique(AbData10$Year)) #make a dataframe to imupt data
names(EvenData10) <- "Year"
EvenData10$Evenness <- 0

for (y in unique(AbData10$Year)){

	AbData10_y <- AbData10[AbData10$Year == y,]		#select data for approprate year
	evennessData_y <- aggregate(AbData10_y$abundance, list(AbData10_y$Species_Identity), sum)      #how many individuals of each species

	H <- diversity(evennessData_y[,2])#Shannon index
	J <- H/log(specnumber(evennessData_y[,2]))#calculate Pileu'e evenness 

	EvenData10$Evenness[EvenData10$Year == y] <- J 

	}


#make code to loop through all data 
#---------------------------------------------------

EvenChange <- list()

i <- 1

for(s in unique(AbData$Study_ID)){


	AbData_s <- AbData[AbData$Study_ID == s,]
	EvenData_s <- data.frame(unique(AbData_s$Year)) #make a dataframe to imupt data
	names(EvenData_s) <- "Year"
	EvenData_s$Evenness <- 0
	EvenData_s$Study_ID <- s

	for (y in unique(AbData_s$Year)){

		AbData_s_y<- AbData_s[AbData_s$Year == y,]		#select data for approprate year
		evennessData_y <- aggregate(AbData_s_y$abundance, list(AbData_s_y$Species_Identity), sum)      #how many individuals of each species

		H <- diversity(evennessData_y[,2])#Shannon index
		J <- H/log(specnumber(evennessData_y[,2]))#calculate Pileu'e evenness 

		EvenData_s$Evenness[EvenData_s$Year == y] <- J 

		}

	EvenChange[[i]] <- EvenData_s

	i <- i +1

}


evennessAbData <- do.call(rbind, EvenChange)
write.csv(evennessAbData, "evennessAbundance.csv")










