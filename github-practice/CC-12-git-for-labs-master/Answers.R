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
LPI.long$country.list <- gsub(",", "", LPI.long$country.list, fixed = TRUE)
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

biome.plots <- LPI.long %>%
  group_by(., genus.species.id, biome) %>% 
  do(mod = lm(scalepop ~ year, data = .)) %>% 
  tidy(mod) %>% 
  select(., genus.species.id, biome, term, estimate) %>% 
  spread(., term, estimate) %>%
  ungroup() %>%
  group_by(., biome) %>%
  do(ggsave(ggplot(.,aes(x = year)) + geom_histogram(colour = "#8B5A00", fill = "#CD8500") + theme_LPI() +
            xlab("Rate of population change (slopes)"), 
            filename = gsub("", "", paste("Biome_LPI/", unique(as.character(.$biome)), ".pdf", sep="")),
            device = "pdf"))

# Challenge 2: Slope estimates and SE vs study duration + marginal histograms ----
# Make a scatterplot of slope estimates and standard errors for all populations vs study duration
# Add histograms of slope estimates and study duration along the margins of the plots

# HINTS:
# You can use the package ggExtra to easily make marginal histograms for ggplot2 objects
(all.slopes <- ggplot(LPI.models, aes(x = lengthyear, y = slope)) +
         geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.4) +
         geom_hline(yintercept = 0, linetype = "dashed") +
         theme_LPI() +
         ylab("Population change\n") +
         xlab("\nDuration (years)"))

all.slopes <- ggExtra::ggMarginal(
  p = all.slopes,
  type = 'histogram',
  margins = 'both',
  size = 5,
  col = 'black',
  fill = 'gray'
)

ggsave(all.slopes, file = "slopes_duration.png", width = 9, height = 7)

# Challenge 3: Histograms for each system type (freshwater, marine, terrestrial) ----
# Make a three panel figure that includes a histogram of slope estimates for each system type

# HINTS:
# The object we made of model outputs (LPI.models) doesn't include the system column
# Bring it back and then you can easily filter by system type and make the histograms

LPI.models2 <- LPI.long %>%
  group_by(., genus.species.id, lengthyear, system, class) %>% 
  do(mod = lm(scalepop ~ year, data = .)) %>%
  mutate(., n = df.residual(mod),
         intercept = summary(mod)$coeff[1],
         slope = summary(mod)$coeff[2],
         intercept_se = summary(mod)$coeff[3],
         slope_se = summary(mod)$coeff[4],
         intercept_p = summary(mod)$coeff[7],
         slope_p = summary(mod)$coeff[8]) %>%
  ungroup() %>%
  mutate(., lengthyear = lengthyear,
         system = system) %>%
  filter(., n > 5) # Remove rows where degrees of freedom <5

# Create objects for each system type
freshwater <- filter(LPI.models2, system == "Freshwater")
marine <- filter(LPI.models2, system == "Marine")
terrestrial <- filter(LPI.models2, system == "Terrestrial")

(hist.fresh <- ggplot(freshwater, aes(x = slope)) + 
    geom_histogram(colour = "#00B2EE", fill = "#00B2EE") + 
    theme_LPI() + 
    labs(x = "Rate of population change (slopes)", title = "(a) Freshwater"))

(hist.mar <- ggplot(marine, aes(x = slope)) + 
    geom_histogram(colour = "#00008B", fill = "#00008B") + 
    theme_LPI() + 
    labs(x = "Rate of population change (slopes)", title = "(b) Marine"))

(hist.terr <- ggplot(terrestrial, aes(x = slope)) + 
    geom_histogram(colour = "#CD950C", fill = "#CD950C") + 
    theme_LPI() + 
    labs(x = "Rate of population change (slopes)", title = "(c) Terrestrial"))

system_panel <- grid.arrange(hist.fresh, hist.mar, hist.terr, ncol = 1)
ggsave(system_panel, file = "system_panel.png", width = 7, height = 15)

# Challenge 4: Population change in mammals, birds, amphibians and reptiles
# Make a four panel figure of scatterplots for the different animal classes
# Slope estimates and SE vs study duration + marginal histograms
# Add animal silhouettes in the top right corner as a legend
# (if you like that sort of stuff)

# HINTS:

# The object we made of model outputs (LPI.models) doesn't include the class column
# Bring it back and then you can easily filter by class and make the scatterplots

# Use the package ggExtra for the marginal histograms

# You can find the animal images in the repo along with the data
# Use the packages png and grid to get the images into a format R recognises
# The functions for that are readPNG() and rasterGrob()

LPI.models2 <- LPI.long %>%
  group_by(., genus.species.id, lengthyear, system, class) %>% 
  do(mod = lm(scalepop ~ year, data = .)) %>%
  mutate(., n = df.residual(mod),
         intercept = summary(mod)$coeff[1],
         slope = summary(mod)$coeff[2],
         intercept_se = summary(mod)$coeff[3],
         slope_se = summary(mod)$coeff[4],
         intercept_p = summary(mod)$coeff[7],
         slope_p = summary(mod)$coeff[8]) %>%
  ungroup() %>%
  mutate(., lengthyear = lengthyear,
         class = class) %>%
  filter(., n > 5) # Remove rows where degrees of freedom <5

# Create objects for each system type
mammals <- filter(LPI.models2, class == "Mammalia")
birds <- filter(LPI.models2, class == "Aves")
amphibians <- filter(LPI.models2, class == "Amphibia")
reptiles <- filter(LPI.models2, class == "Reptilia")

# Load packages for adding images
packs <- c("png","grid")
lapply(packs, require, character.only = TRUE) 

icon_m <- readPNG("tiger.png") 
icon_m <- rasterGrob(icon_m, interpolate = TRUE)

(slopes.mammals <- ggplot(mammals, aes(x = lengthyear, y = slope)) +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("\nPopulation change\n") +
    xlab("\nDuration (years)\n") + labs(title = "\n(a) Mammals\n") +
    annotation_custom(icon_m, xmin = 30, xmax = 40, ymin = 0.06, ymax = 0.16))

slopes.mammals <- ggExtra::ggMarginal(
  p = slopes.mammals,
  type = 'histogram',
  margins = 'y',
  size = 5,
  col = 'black',
  fill = 'gray'
)

icon_b <- readPNG("bird.png") 
icon_b <- rasterGrob(icon_b, interpolate = TRUE)

(slopes.birds <- ggplot(birds, aes(x = lengthyear, y = slope)) +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("\nPopulation change\n") +
    xlab("\nDuration (years)\n") + labs(title = "\n(b) Birds\n") +
    annotation_custom(icon_b, xmin=34, xmax=44, ymin=0.1, ymax=0.18))

slopes.birds <- ggExtra::ggMarginal(
  p = slopes.birds,
  type = 'histogram',
  margins = 'y',
  size = 5,
  col = 'black',
  fill = 'gray'
)

icon_a <- readPNG("frog.png") 
icon_a <- rasterGrob(icon_a, interpolate = TRUE)

(slopes.amphibians <- ggplot(amphibians, aes(x = lengthyear, y = slope)) +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("\nPopulation change\n") +
    xlab("\nDuration (years)\n") + labs(title = "\n(c) Amphibians\n") +
    annotation_custom(icon_a, xmin=24, xmax=28, ymin=0.06, ymax=0.16))

slopes.amphibians <- ggExtra::ggMarginal(
  p = slopes.amphibians,
  type = 'histogram',
  margins = 'y',
  size = 5,
  col = 'black',
  fill = 'gray'
)

icon_r <- readPNG("snake.png") 
icon_r <- rasterGrob(icon_r, interpolate=TRUE)

(slopes.reptiles <- ggplot(reptiles, aes(x = lengthyear, y = slope)) +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se), alpha = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("\nPopulation change\n") +
    xlab("\nDuration (years)\n") + labs(title = "\n(d) Reptiles\n") +
    annotation_custom(icon_r, xmin=30, xmax=38, ymin=0.08, ymax=0.16))

slopes.reptiles <- ggExtra::ggMarginal(
  p = slopes.reptiles,
  type = 'histogram',
  margins = 'y',
  size = 5,
  col = 'black',
  fill = 'gray'
)
class_panel <- grid.arrange(slopes.mammals, slopes.birds, slopes.amphibians, slopes.reptiles, ncol = 2)
ggsave(class_panel, file="class_panel.png", width = 12, height = 10)

# Bonus graphs
(slopes_freshwater <- ggplot(freshwater, aes(x = lengthyear, y = slope)) +
    geom_point() +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("Population change\n") +
    xlab("\nDuration (years)") + labs(title = "Freshwater"))

marine <- filter(LPI_models2, system == "Marine")
(slopes_marine <- ggplot(marine, aes(x = lengthyear, y = slope)) +
    geom_point() +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("Population change\n") +
    xlab("\nDuration (years)") + labs(title = "Marine"))

terrestrial <- filter(LPI_models2, system == "Terrestrial")
(slopes_terrestrial <- ggplot(terrestrial, aes(x = lengthyear, y = slope)) +
    geom_point() +
    geom_pointrange(aes(ymin = slope - slope_se, ymax = slope + slope_se)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_LPI() +
    ylab("Population change\n") +
    xlab("\nDuration (years)") + labs(title = "Terrestrial"))

system_panel2 <- grid.arrange(slopes_freshwater, slopes_marine, slopes_terrestrial, ncol = 1)
ggsave(system_panel2, file="system_panel2.png", width = 7, height = 15)
