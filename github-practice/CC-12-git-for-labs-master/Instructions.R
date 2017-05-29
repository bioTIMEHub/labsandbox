# Practice script for GitHub lab workshop
# Gergana Daskalova & Isla Myers-Smith
# 29-05-2017

# Packages ----
library(readr)
library(tidyr)
library(dplyr)
library(broom)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(gridExtra)

# Load data ----
load("LPIdata_Feb2016.RData")
View(head(LPIdata_Feb2016))

# Format data for analysis ----

# Transform from wide to long format
LPI.long <- gather(data = LPIdata_Feb2016, key = "year", value = "pop", select = 26:70)

# Get rid of the X in front of years
LPI.long$year <- parse_number(LPI.long$year)

# Rename variable names for consistency
names(LPI.long)
names(LPI.long) <- tolower(names(LPI.long))
names(LPI.long)

# Create new column with genus, species and id to differentiate between different populations
LPI.long$genus.species.id <- paste(LPI.long$genus, LPI.long$species, LPI.long$id, sep = "_")

# Check data are displayed fine
View(LPI.long[c(1:5,500:505,1000:1005),])  
# You can use [] to subset data frames [rows, columns]
# If you want all rows/columns, add a comma in the row/column location

# Get rid of strange characters like " / "
LPI.long$country_list <- gsub(",", "", LPI.long$country.list, fixed = TRUE)
LPI.long$biome <- gsub("/", "", LPI.long$biome, fixed = TRUE)

# Examine the tidy data frame
View(head(LPI.long))

# Data manipulation ----

# Remove duplicate rows
LPI.long <- distinct(LPI.long)

# Remove missing / infinite data
LPI.long <- filter(LPI.long, is.finite(pop))

# Keep species with >5 years worth of data 
# Calculate length of monitoring and scale population trend data
LPI.long <- LPI.long %>%
  group_by(., genus.species.id) %>%  # group rows so that each group is one population
  mutate(., maxyear = max(year), minyear = min(year)) %>%  # Create columns for the first and most recent years that data was collected
  mutate(., lengthyear = maxyear-minyear) %>%  # Create a column for the length of time data available
  mutate(., scalepop = (pop-min(pop))/(max(pop)-min(pop))) %>%
  filter(., is.finite(scalepop)) %>%
  filter(., lengthyear > 5) %>%
  ungroup(.)

# Modelling population change over time ----
# Run linear models of abundance trends over time for each population and extract model coefficients
LPI.models <- LPI.long %>%
  group_by(., genus.species.id, lengthyear) %>% 
  do(mod = lm(scalepop ~ year, data = .)) %>%  # Create a linear model for each group
  mutate(., n = df.residual(mod),  # Create columns: degrees of freedom
         intercept = summary(mod)$coeff[1],  # intercept coefficient
         slope = summary(mod)$coeff[2],  # slope coefficient
         intercept_se = summary(mod)$coeff[3],  # standard error of intercept
         slope_se = summary(mod)$coeff[4],  # standard error of slope
         intercept_p = summary(mod)$coeff[7],  # p value of intercept
         slope_p = summary(mod)$coeff[8]) %>%  # p value of slope
  ungroup() %>%
  mutate(., lengthyear = lengthyear) %>%
  filter(., n > 5) # Remove rows where degrees of freedom <5

# Visualising model outputs ----

# Setting a custom ggplot2 function
theme_LPI <- function(){
  theme_bw()+
    theme(axis.text.x=element_text(size=12, angle=45, vjust=1, hjust=1),
          axis.text.y=element_text(size=12),
          axis.title.x=element_text(size=14, face="plain"),             
          axis.title.y=element_text(size=14, face="plain"),             
          panel.grid.major.x=element_blank(),                                          
          panel.grid.minor.x=element_blank(),
          panel.grid.minor.y=element_blank(),
          panel.grid.major.y=element_blank(),  
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          plot.title = element_text(size=20, vjust=1, hjust=0.5),
          legend.text = element_text(size=12, face="italic"),          
          legend.title = element_blank(),                              
          legend.position=c(0.9, 0.9))
}

# Challenge 1: Slope estimates for each biome ----
# Make histograms of slope estimates for each biome and save them in a folder called Biome_LPI

# HINTS:
# You need to make the folder where you save the files first
# You can use dplyr to group by biome and run simple linear models (e.g. population size through time)
# With broom you can get a dataframe with model outputs
# Use do() before a function call within a pipe to make the function run for all groupings in the pipe

# Challenge 2: Slope estimates and SE vs study duration + marginal histograms ----
# Make a scatterplot of slope estimates and standard errors for all populations vs study duration
# Add histograms of slope estimates and study duration along the margins of the plots
# Save your graph in the github-practice folder

# HINTS:
# You can use the package ggExtra to easily make marginal histograms for ggplot2 objects

# Challenge 3: Histograms for each system type (freshwater, marine, terrestrial) ----
# Make a three panel figure that includes a histogram of slope estimates for each system type
# Save your figure in the github-practice folder

# HINTS:
# The object we made of model outputs (LPI.models) doesn't include the system column
# Bring it back and then you can easily filter by system type and make the histograms

# Challenge 4: Population change in mammals, birds, amphibians and reptiles ----
# Make a four panel figure of scatterplots for the different animal classes
# Slope estimates and SE vs study duration + marginal histograms
# Add animal silhouettes in the top right corner as a legend
# (if you like that sort of stuff)
# Save your panel in the github-practice folder

# HINTS:

# The object we made of model outputs (LPI.models) doesn't include the class column
# Bring it back and then you can easily filter by class and make the scatterplots

# Use the package ggExtra for the marginal histograms

# You can find the animal images in the repo along with the data
# Use the packages png and grid to get the images into a format R recognises
# The functions for that are readPNG() and rasterGrob()