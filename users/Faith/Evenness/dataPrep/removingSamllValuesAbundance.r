library(parallel)
library(foreach)
library(doParallel)

setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
DataAbundanceOnly <- read.csv("DataAbundanceBeforeRemoval.csv")

data278 <- DataAbundanceOnly[DataAbundanceOnly$STUDY_ID == 278,]
summary(data278)
head(data278)

head(DataAbundanceOnly)
LargeSamples<-list() 

no_cores <- detectCores()-1
cl <- makeCluster(no_cores)

cl<-makeCluster(detectCores()-1)
registerDoParallel(cl)

unique(DataAbundanceOnly[DataAbundanceOnly$STUDY_ID == 284, 3])

LargeSamples <- foreach(b = unique(DataAbundanceOnly$STUDY_ID)) %dopar% {

	datab <- DataAbundanceOnly[DataAbundanceOnly$STUDY_ID == b,]

	nsampYear  <- data.frame(matrix(NA, nrow = length(unique(datab$YEAR)), ncol = 2))
	nsampYear[,1] <- unique(datab$YEAR)

		for (year in unique(datab$YEAR)) {
			nsamp <-  length(unique(datab$SAMPLE_DESC[datab$YEAR==year]))
			nsampYear[,2][nsampYear[,1] == year] <- nsamp
			colnames(nsampYear) <- c("Year", "Sample_length") 
		}

	smallSample <- mean(nsampYear$Sample_length)/2 # calculating half of the average number of samples 

	#remove years from datab whos number of samples is less than half average (smallSample)

	yearsNotSmall <- subset(nsampYear, nsampYear$Sample_length > smallSample)

	databNoSmall <- datab[datab$YEAR %in% yearsNotSmall$Year, ]

}
stopCluster(cl)
dataAb <-  do.call("rbind", LargeSamples)

head(dataAb)

write.csv(dataAb, "AbundanceDataNoSmallVales.csv")



