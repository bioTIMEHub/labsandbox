setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\data")
library(ggplot2)
library(ggExtra)
library(ellipse)

evennessSim <- read.csv("evennessChangeSimp.csv")
evenness <- read.csv("evennessChange.csv")
Abundance <- read.csv("AbundanceCapacityChange.csv")
AbDominance <- read.csv("ActualAbundanceChange.csv")
RelDominance <- read.csv("PercentAbundanceChange.csv")
rareity <- read.csv("rarityChange.csv")
rareity2 <- read.csv("rarityChangeSingelton.csv")
SRChange <- read.csv("SpeciesRichnessChange.csv")

##renaming columns so they make sense 
names(Abundance)[4] <- "AbundanceAll"
names(evennessSim)[4] <- "EvennessSim1"
names(evenness)[4] <- "Evenness"
names(AbDominance)[4] <- "AbundanceDom"
names(RelDominance)[4] <- "relativeDom"
names(rareity)[4] <- "Rareity"
names(rareity2)[4] <- "Rareity2"

##merging all data
AbEv <- merge(Abundance, evenness, by = "Study")
AbRel <- merge(AbDominance, RelDominance, by = "Study")
RareChange <- merge(rareity,AbRel, by = "Study")
RareChange2 <- merge(rareity2,RareChange, by = "Study")
Allchange1 <- merge(AbEv,RareChange2, by = "Study")
AllSR <- merge(SRChange,Allchange1, by.x = "Study_ID", by.y = "Study")
Allchange <- merge(evennessSim, AllSR,  by.y = "Study_ID", by.x = "Study")
names(Allchange)

#selecting only columns with coefficients of change, and not all teh
#extra info like taxa or ;ength of study
selectData <- Allchange[,c(1,4, 13, 18, 29, 38, 44, 50, 60, 5, 61, 62, 41)]
head(selectData)

##plot scatter plot with marginal histograms
#evennes-assemblage size
evAss <- ggplot(selectData, aes( AbundanceAll, Evenness)) + 
	geom_point() + 
	theme_classic()
ggMarginal(evAss, type = "histogram")

#chaneg taxa to just invert/vert/plant
levels(selectData$TAXA.x)
levels(selectData$TAXA.x)[1] <- "Multipe"
levels(selectData$TAXA.x)[2] <- "Vertebrate"
levels(selectData$TAXA.x)[4] <- "Vertebrate"
levels(selectData$TAXA.x)[4] <- "Vertebrate"
levels(selectData$TAXA.x)[5] <- "Vertebrate"
levels(selectData$TAXA.x)[6] <- "Vertebrate"
levels(selectData$TAXA.x)[3] <- "Multipe"
selectData <- unique(selectData)

#check soem stats on the data
#------------------------------

aggdata <- aggregate(selectData$Study, by=list(selectData$REALM, selectData$PROTECTED_AREA.y.1, selectData$TAXA.x),
  FUN=length)
names(aggdata) <- c("Realm", "Preotected_area", "Taxa", "Count")

#calculating correlaion matrix
#-------------------------------
head(selectData)

selectDataCor <- selectData[,c(2:4)]
head(selectDataCor)
row.names(selectDataCor) <- selectData$Study

corMat <- cor(selectDataCor)

plorEllipse <- data.frame(ellipse(corMat[,c(1:2)]))




#EVENESS and ABUNDANCE 
#----------------------

#evenness and abndance
evAss <- ggplot(selectData, aes(AbundanceAll, Evenness)) + 
	geom_point() + 
	theme_classic()+ 
	scale_colour_manual(values = c("DarkViolet", "Orange","LightBlue", "Green"))+
	labs(x = "Assembalge Size Change", y = "Evenness Change")


ggMarginal(evAss, type = "histogram")
	

ggMarginal(evAssT, type = "histogram")
#evennes-assemblage size TAXA
evAssT <- ggplot(selectData, aes(AbundanceAll, Evenness, col = TAXA.x)) + 
	geom_point() + 
	theme_classic()+ 
	scale_colour_manual(values = c("DarkViolet", "Orange","LightBlue", "Green"))+
	labs(x = "Assembalge Size Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)

ggMarginal(evAssT, type = "histogram")


#evennes-assemblage size Protected Area
evAssPA <- ggplot(selectData, aes( AbundanceAll, Evenness, col = PROTECTED_AREA.y.1)) + 
	geom_point() + 
	theme_classic() + 
	labs(x = "Assembalge Size Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAssPA, type = "histogram")

#evennes-assemblage size LengthTime
evAssL <- ggplot(selectData, aes( AbundanceAll, Evenness, size = LengthYears.y)) + 
	geom_point() + 
	theme_classic()+ 
	labs(x = "Assembalge Size Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAssL, type = "histogram")

#evennes-assemblage size REALM
evAssR <- ggplot(selectData, aes( AbundanceAll, Evenness, col = REALM)) + 
	geom_point() + 
	theme_classic()+ 
	scale_colour_manual(values = c("Orange","LightBlue", "Green"))+
	labs(x = "Assembalge Size Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAssR, type = "histogram")

#EVENNESS and SPECIES RICHNESS
#-----------------------------

evAb <- ggplot(selectData, aes( SRslope, Evenness)) + 
	geom_point() + 
	theme_classic()+ 
	labs(x = "Species Richness Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAb, type = "histogram")

#evennes-species richness TAXA
evAbT <- ggplot(selectData, aes( SRslope, Evenness, col = TAXA.x)) + 
	geom_point() + 
	theme_classic()+ 
	scale_colour_manual(values = c("DarkViolet", "Orange","LightBlue", "Green"))+
	labs(x = "Species Richness Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAbT, type = "histogram")

#evennes-species richness Protected Area
evAbPA <- ggplot(selectData, aes( SRslope, Evenness, col = PROTECTED_AREA.y.1)) + 
	geom_point() + 
	theme_classic()+ 
	labs(x = "Species Richness Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAbPA, type = "histogram")

#evennes-assemblage size REALM
evAbSR <- ggplot(selectData, aes(SRslope, Evenness, col = REALM)) + 
	geom_point() + 
	theme_classic()+ 
	scale_colour_manual(values = c("Orange","LightBlue", "Green"))+
	labs(x = "Species Richness Change", y = "Evenness Change")+
	stat_ellipse(type = "norm", level = 0.95)
ggMarginal(evAbSR, type = "histogram")

###corrilation between slopes
selectData$SRslope
cor(selectData$SRslope, selectData$Evenness,  method = "spearman")
cor(selectData$AbundanceAll, selectData$Evenness,  method = "pearson")
cor(selectData$AbundanceAll, selectData$SRslope,  method = "pearson")

names(selectData)
head(selectData)

##PCR fisrt attempt
library(ggfortify)
library(cluster)
library(pvclust)

autoplot(prcomp(scale(selectData[,c(2, 3, 4, 6, 8, 9)])), data = selectData, colour = 'PROTECTED_AREA.y.1', loadings = TRUE,
	loadings.label = TRUE)

autoplot(prcomp(scale(selectData[,c(2, 3, 4, 12)])), data = selectData, colour = 'PROTECTED_AREA.y.1', loadings = TRUE,
	loadings.label = TRUE)

plot(selectData[,c(2, 3, 4, 6, 12)])

ScaleData <- data.frame(scale(selectData[,c(2, 3, 4, 6)]))
plot(ScaleData)

autoplot(pam(scale(selectData[,c(2, 3, 4, 6, 12)]), 3), frame = TRUE, frame.type = 'norm')

clusterp <- hclust(dist(selectData[,c(2, 3, 4, 6)]))
plot(clusterp)
#dendrogram 
mydata <- scale(selectData[,c(2, 3, 4, 12)])

#determine number of clusters 
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(mydata,
   centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters",
  ylab="Within groups sum of squares") 

# K-Means Cluster Analysis
fit <- kmeans(mydata, 3) # 3 cluster solution
# get cluster means
aggregate(mydata,by=list(fit$cluster),FUN=mean)
# append cluster assignment
mydata <- data.frame(mydata, fit$cluster) 

#plot clusters on a scatter graph
denPlot <- ggplot(data = mydata, aes(x = SRslope, y = EvennessSim1, col = as.factor(fit.cluster)))
denPlot + geom_point()

# Ward Hierarchical Clustering
d <- dist(mydata, method = "euclidean") # distance matrix
fit <- hclust(d, method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=3) # cut tree into 3 clusters
# draw dendogram with red borders around the 3 clusters
rect.hclust(fit, k=3, border="red") 

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
fit <- pvclust(mydata, method.hclust="ward",
   method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95) 


head(mydata)
#do any factors explain clustering?
mydata <- data.frame(mydata, selectData$REALM, selectData$TAXA.x)
names(mydata)[6] <- "REALM"
names(mydata)[7] <- "TAXA"
plot(mydata$REALM, mydata$fit.cluster)

tmydata <- t(mydata[,1:4])

library(pvclust)
fit <- pvclust(tmydata, method.hclust="ward",
   method.dist="euclidean")
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
pvrect(fit, alpha=.95)



# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 3)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster)
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,
   labels=2, lines=0)

# Centroid Plot against 1st 2 discriminant functions
library(fpc)
plotcluster(mydata, fit$cluster) 


####################
## Solving the question:

# loading the package
library(dendextend)
# Assigning the labels of dendrogram object with new colors:
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
# Plotting the new dendrogram
plot(dend)




