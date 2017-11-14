#Function to rarefy data once and report individual abundances per species 
#-------------------------------------------

rarefyOnce <- function(Year, SampleID, Biomass, Species){

  #######################################################################
  # takes as input a  Year, SampleID, Biomass, and Species 
  # which should be in dataframe so that elements match
  ###########################################################################

	#select year
	years <- unique(Year)#which years are in this data set
	DataStudyYear <- Year[Year == years[1]]

	#getting vector with number of samples per year
	nsamples<-c()
 	for(y in unique(Year)){
   		nsamples<-c(nsamples, length(unique(SampleID[Year==y])))
  	}
	minsample<-min(nsamples)
	t<-1
	BiomassSpecies <- list()
	for(y in unique(Year)){	
		rareBiomass <- c()
		samps<-unique(SampleID[Year==y])
      	# re-sampling samples to equalize number of samples
    		sam<-as.character(sample(samps,minsample, replace=T))
		# Get species capacity and identity from that year
		SampleID <- as.character(SampleID)
		rareBio <- Biomass[SampleID %in% sam]
		rareSpecies <- Species[which (SampleID %in% sam)]
		bioSpecies <- data.frame(cbind(rareBio,y))
		bioSpecies$Spcies <- rareSpecies
		bioSpecies <- bioSpecies[complete.cases(bioSpecies),]
		BiomassSpecies[[t]] <- bioSpecies
		t<-t+1
	}
	RareBio <- do.call(rbind, BiomassSpecies)#bind all lists together
  	return(RareBio)
   }

