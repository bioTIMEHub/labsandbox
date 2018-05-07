setwd("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\NullModel_correctsubset")

#code to compare the Z scores of each of dominance change and assemblage chaneg measurements together.
#data was writen in the "AssemblageSize.r, ActualDominance.r and PercentDominance.r sripts. Each number
#is (slope - meanNullModelSlope)/NullModelSlope_sandardDeviation

library(ggplot2)
library(dplyr)
library(cowplot)
library(rgl)#3d plot
library(car)#3d plot
library(grid)
library(gridExtra)

#i used teh RS model from teh previous null modle because teh more recent one failed to converge. This
#is just for plotting purposes. 
SRZscore <- read.csv("C:\\Users\\faj\\Documents\\OneDrive for Business\\research\\ch3 - metaanalysis\\Oct2017Analysis\\Null Model\\SpecRich_Zscores.csv")
AssSizeZscore <- read.csv("AssSize_Zscores.csv")
RelDomZscore <- read.csv("RelDom_Zscores.csv")
ActDomZscore <- read.csv("ActDom_Zscores.csv")

SRZscore$X <- NULL
AssSizeZscore$X <- NULL
RelDomZscore$X <- NULL
ActDomZscore$X <- NULL
names(RelDomZscore)[2] <- "relDomZscores"
nrow(SRZscore)


SRass<- merge(SRZscore, AssSizeZscore, by = "Study_ID")
SRassRel<- merge(RelDomZscore, SRass, by = "Study_ID")
AllZscores<- merge(SRassRel, ActDomZscore, by = "Study_ID")
head(AllZscores)

#remove study 360 because although it is listed as having count data, all 
#abundances are 0. Thsi corrected in teh new query
#We also remove study 248 because it is too experimental 
#merge info together

AllZscores[AllZscores$Study_ID %in% c("360", "248"),]


#plot each variabel against year other

SRRel <- ggplot(AllZscores,aes(x = SRZscores, y = relDomZscores)) +
	geom_point()+
	theme_classic()+
	labs(x = "", y = "Relative Dominance Z score")
SRAbs <- ggplot(AllZscores,aes(x = SRZscores, y = ActDomZscores)) +
	geom_point()+
	theme_classic()+
	labs(x = "Species Richness Z score", y = "Absolute Dominance Z score")

AssRel <- ggplot(AllZscores,aes(x = AssSizeZscores, y = relDomZscores)) +
	geom_point()+
	theme_classic()+
	labs(x = "", y = "")
cor(AllZscores$AssSizeZscores, AllZscores$relDomZscores, method = "pearson")
cor(AllZscores$AssSizeZscores, AllZscores$relDomZscores, method = "spearman")

AssAbs <- ggplot(AllZscores,aes(x = AssSizeZscores, y = ActDomZscores)) +
	geom_point()+
	theme_classic()+
	labs(x = "Assemblage Size Z score", y = "")
cor(AllZscores$AssSizeZscores, AllZscores$ActDomZscores, method = "pearson")
cor(AllZscores$AssSizeZscores, AllZscores$ActDomZscores, method = "spearman")

plot_grid(SRRel, AssRel, SRAbs, AssAbs, labels = c("A", "B", "C", "D"))

#3d plot of z scores 
#--------------
scatter3d(ActDomZscores ~ AssSizeZscores + SRZscores, data=AllZscores, fit="smooth", 
	xlab = "Assemblage size Change", ylab = "Absolute dominance change", zlab = "Species richness Change")
scatter3d(relDomZscores ~ AssSizeZscores + SRZscores, data=AllZscores, fit="smooth", 
	xlab = "Assemblage size Change", ylab = "relative dominance change", zlab = "Species richness Change")


#making to scatter plots for a presentation

AssRel2 <- ggplot(AllZscores,aes(x = AssSizeZscores, y = relDomZscores)) +
	geom_point(col = "red")+
	theme_classic()+
	labs(x = "Assemblage Size Z score", y = "Relative dominance z score")+
	theme(text = element_text(size=20))+
  	geom_text(aes(x = -8, y = 110), label = "B)", size = 6)

AssAbs2 <- ggplot(AllZscores,aes(x = AssSizeZscores, y = ActDomZscores)) +
	geom_point(col = "blue")+
	theme_classic()+
	labs(x = "", y = "Absolute dominance z score")+
	theme(text = element_text(size=20))+
  	geom_text(aes(x = -8, y = 9), label = "A)", size = 6)


bothPlot <- grid.arrange(AssAbs2, AssRel2, ncol = 1)

#ggsave("Zscores.jpeg", plot = bothPlot, device = "jpeg", width = 15, height = 22, units = c("cm"))



