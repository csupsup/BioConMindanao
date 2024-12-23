---
title: "Supplementary FS3"
date: "16 December 2024"
output:
  html_document:
    df_print: paged
  pdf_document: default
indent: no
---
```{r setup, include=FALSE, message=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
```
This file documents our reanalysis of the dataset used to examine biodiversity and conflict relationships in the\ 
Southern Philippines presented in this paper https://doi.org/10.1038/s44185-024-00044-8.\ 
Here, we consider spatial autocorrelation.

**Files needed**

```{r files, eval = FALSE}
##1. "cdh" folder containing the raster files of tree cover, tree density and forest canopy height
##2. "landcover" folder containing the raster files of 1988 and 2020 land cover types
##3. "mobios_tab.csv" or the biodiversity data
##4. "conflict_tab_89-21.csv" or the conflict data
##5. "functions.R" for calculating non-significant autocorrelated average distance
##6. "1988_lc_classmap.csv" 1988 land cover types code
##7. "2020_lc_classmap.csv" 2020 land cover types code
##8. "mindanao_adm_geo.gpkg" Mindanao boundary
```

**Load required packages**

```{r load package, message=FALSE}
## make sure to install all packages including required dependencies prior to use
re.lib <- c("dplyr", "tidyverse", "ggplot2", "patchwork", "GeoThinneR", "ecodist", 
                        "sf", "sp", "spdep", "raster", "MASS", "terra")

new_packages <- re.lib[!(re.lib %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)

# load packages
invisible(lapply(re.lib, library, character.only = TRUE))
```

**Prepare datasets**

```{r read data}
## biodiversity data obtained from https://doi.org/10.15468/rtedgk
## conflict data obtained from https://data.humdata.org/dataset/ucdp-data-for-philippines?
## note that only a subset of this dataset (i.e., those from Mindanao and adjacent islands) was used

## read biodiversity data
bio <- read.csv("data/mobios_tab.csv", header = TRUE, sep = ",")
#subset data
bio.sub <- subset(bio, select = c(12, 23, 20, 15, 16))
bio.sub <- bio.sub %>% filter(county!="") #remove rows with no county information
bio.sub <- bio.sub %>% filter(class!="Malacostraca") #exclude this taxon
bio.sub <- bio.sub %>% filter(class!="Bivalvia") #exclude this taxon
bio.sub <- bio.sub %>% filter(class!="Gastropoda")#exclude this taxon

## check province name
## there are rows where the name of the province is uncertain as indicated by "/"
## "Zamboanga¬†del Norte" was changed to "Zamboanga del Norte" prior to import of data
unique(bio.sub$county)

## remove uncertain province names
bio.sub <- bio.sub %>% filter(!grepl('/', county))
head(bio.sub)

## write data
write.csv(bio.sub, "bio.sub.csv")

## read conflict data
con <- read.csv("data/conflict_tab_89-21.csv", header = TRUE, sep = ",")

## subset to get revelant data
con.sub <- subset(con, select = c(1,3,15,18,30,32,33)) 
head(con.sub)

## export filtered conflic data
write.csv(con.sub, "con.sub.csv")
```

**Address spatial autocorrelation in both biodiversity and conflict data**

```{r spatialautocorrelation}
## this part was not done in the paper, which may affect the analysis
## first, filter species with duplicate coordinates
bio.unique <- bio.sub %>%
  group_by(scientificName) %>%
  distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

## then find the non-significant autocorrelated average distance to thin the data\ 
## based on the relationship between land cover types

## use the spatialautocorrelation function 
## from https://github.com/jorgeassis/spatialAutocorrelation
source("data/functions.R")

## read data
files <- list.files("data/landcover/", pattern = "tif", full.names=TRUE)
raster.files <- lapply(files, raster)

raster.files

## make rasters to have common extent
common_extent <- do.call(union, lapply(raster.files, extent))

raster.files <- lapply(raster.files, function(r) {
    extend(r, common_extent)
    extent(r) <- common_extent 
    r
})

## check extent
lapply(raster.files, extent)

## resample raster files to have a uniform dimension/resolution
standard <- raster.files[[1]]
raster.files <- lapply(raster.files, function(r) {
    resample(r, standard, method = "bilinear")
})

## check resolution
lapply(raster.files, res)

## stack the raster files
lc.data <- stack(raster.files)

## get longitude and latitude columns
bio.coor <- as.data.frame(subset(bio.unique, select = c("decimalLongitude", "decimalLatitude")))

## define the distance class (in km)
autocorrelationClassDistance <- 2

## define the maximum distance (in km)
autocorrelationMaxDistance <- 10

## define the significance level of the test
autocorrelationSignif <- 0.05

distanceUncorr <- data.frame(Predictor=names(lc.data),Distance=NA)

for( i in 1:length(names(lc.data)))
  distanceUncorr[i,2] <- spatialAutocorrelation(occurrenceRecords=bio.coor,subset(lc.data,i), 
    autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif)

meanCorrDistance <- mean(distanceUncorr[,2])

meanCorrDistance

## thin every species data using the mean distance
## species name
sp.name <- unique(bio.unique$scientificName)

## create an empty list to store the results
thin.data <- list()

## go over each species
for (spp in sp.name) {
  bio.spp <- subset(bio.unique, scientificName == spp)
  bio.thin <- thin_points(
    data = bio.spp,
    long_col = "decimalLongitude", 
    lat_col = "decimalLatitude",
    method = "brute_force",
    thin_dist = meanCorrDistance,  # thinning distance in km
    trials = 1, # number of replicates
    all_trials = TRUE, # return all trials
    seed = 123 # seed for reproducibility
  )
  
  ## convert results to a data frame
  bio.thin.df <- as.data.frame(bio.thin[1])
  
  ## add the species name
  bio.thin.df$scientificName <- spp
  
  ## append the data to the list
  thin.data[[spp]] <- bio.thin.df
}

## combine all thinned
bio.thin.df <- do.call(rbind, thin.data)

## write to a csv
write.csv(bio.thin.df, "bio.thin.csv")

## convert thin data to sf object
bio.thin.sf <- st_as_sf(bio.thin.df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

## plot biodiversity data
## read Mindanao admin boundary
ph.poly <- st_read("data/mindanao_adm_geo.gpkg") 

## plot 
ggplot() +
  geom_sf(data = ph.poly) +
  geom_sf(data = bio.thin.sf, color = "blue", size = 2) +
  theme_minimal() +
  labs(title = "Thinned Biodiversity Data")

## check the number of records (original vs thinned data)
nrow(bio.sub) # original data
nrow(bio.thin.df) # thinned data

## apply the same apporach to conflict data
con.unique <- con.sub %>%
  group_by(PROVINCE) %>%
  distinct(longitude, latitude, .keep_all = TRUE)

## get longitude and latitude columns
con.coor <- as.data.frame(subset(con.unique, select = c("longitude", "latitude")))

## define the distance class
autocorrelationClassDistance <- 2

## define the maximum distance (in km)
autocorrelationMaxDistance <- 10

## define the significance level of the test
autocorrelationSignif <- 0.05

distanceUncorr <- data.frame(Predictor=names(lc.data),Distance=NA)

for( i in 1:length(names(lc.data)))
  distanceUncorr[i,2] <- spatialAutocorrelation(occurrenceRecords=con.coor,subset(lc.data,i), 
    autocorrelationClassDistance,autocorrelationMaxDistance,autocorrelationSignif)

meanCorrDistance <- mean(distanceUncorr[,2])

meanCorrDistance

## thin every province data using the mean distance
## use province name
prov.name <- unique(con.unique$PROVINCE)

## create an empty list to store the results
thin.data <- list()

## go over each province
for (prov in prov.name) {
  con.prov <- subset(con.unique, PROVINCE == prov)
  con.thin <- thin_points(
    data = con.prov,
    long_col = "longitude", 
    lat_col = "latitude",
    method = "brute_force",
    thin_dist = meanCorrDistance,  # thinning distance in km
    trials = 1, # number of replicates
    all_trials = TRUE, # return all trials
    seed = 123 # seed for reproducibility
  )
  
  ## convert results to a data frame
  con.thin.df <- as.data.frame(con.thin[1])
  
  ## add the province name
  con.thin.df$PROVINCE <- prov
  
  ## append the data to the list
  thin.data[[prov]] <- con.thin.df
}

## combine all thinned
con.thin.df <- do.call(rbind, thin.data)

## write to a csv
write.csv(con.thin.df, "con.thin.csv")

## convert thin data to sf object
con.thin.sf <- st_as_sf(con.thin.df, coords = c("longitude", "latitude"), crs = 4326)

## plot conflict data
## read Mindanao admin boundary
ph.poly <- st_read("data/mindanao_adm_geo.gpkg") 

## plot 
ggplot() +
  geom_sf(data = ph.poly) +
  geom_sf(data = con.thin.sf, color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Thinned Conflict Data")

## check the number of records (original vs thinned data)
nrow(con.sub) # original data
nrow(con.thin.df) # thinned data
```

**Count recorded species per taxon in each point and within province**

```{r count species, out.width="70%"}
## the procedures for conducting species counts in each site were not mentioned in the paper
## according to the paper, analysis was done at the provincial level
## there could be two ways that species counts were done (see image below)
## this part is extremely important, which the authors failed to discuss
## the left panel of the figure shows that species counts were aggregated at the provincial level 
## each point (within the province) receives the same value of species count (or species richness)
## the right panel shows that each point has a unique number of species count
## note that these values are crucial for other downstream analyses 
## particularly when relationships of species counts 
## and distance to the nearest conflict sites are considered
## we generated the values for these two scenarios below 
## and test whether any of them would match the results presented in the paper

knitr::include_graphics("images/bioconflict.png")

## count unique records (i.e., species richness) per point
bio.data.sum <- bio.thin.df %>%
              group_by(county, class, decimalLatitude, decimalLongitude) %>%
              mutate(NumbSp = n_distinct(scientificName))

head(bio.data.sum)

## write to a csv
write.csv(bio.data.sum, "bio.data.sum.csv")

## count species records per province (this lumps all the data within each province, 
## thus each point will have a single species richness)
prov.bio.data.sum <- bio.thin.df %>%
                          group_by(county, class) %>%
                          mutate(NumbSp = n_distinct(scientificName))
head(prov.bio.data.sum)

## write to a csv
write.csv(prov.bio.data.sum, "prov.data.sum.csv")
```

**Calculate distance of species record to the nearest conflict sites using QGIS**

```{r cal distance,  out.width="50%"}
## in QGIS, import the "bio.data.sum.csv" and "con.thin.csv" from R as delimited text files. 
## export these files as GeoPackage or Shapefile to convert the data to meters prior to calculation
## use the CRS UTM Zone 51N
knitr::include_graphics("images/qgis_geo.png")
knitr::include_graphics("images/qgis_utm.png")
## import back the newly saved data to QGIS
## get the distance of nearest conflict areas using the "Joint attributes by nearest" plugin
## QGIS > Processing Toolbox > Join attributes by nearest
knitr::include_graphics("images/qgis_join.png")
## distance can be calculated as well without joining the attributes using the "Distance to nearest hub"
## after calculating the distance, a new layer will be created. 
## the new layer contains the calculated distance in the attribute table
## this also includes the x and y coordinates of the nearest conflict site
## then export this new layer as a csv file so we can use the data in R for other analysis
## repeat the entire process for provincial data (i.e., prov.bio.data.sum.csv) if necessary
## figure below shows the nearest conflict point to species record (black line)
## notice the remainder of conflict points are not included
knitr::include_graphics("images/hub.png")
```

**Calculate average distance**

```{r cal mean distance}
## from now on, we'll solely use specie richness data per point
## read the file with calculated distance from QGIS
bio.con.dist <- read.csv("data/bio.con.dist.csv", header = TRUE, sep = ",")

## get average distance by province for each taxon (i.e., class column)
mean.dist <- bio.con.dist %>%
                group_by(county, class, NumbSp, decimalLatitude, decimalLongitude) %>%
                summarize(meanDistance = mean(distance, na.rm = TRUE))

head(mean.dist)

## write to a csv
write.csv(mean.dist, "mean.distance.csv")
```

**Plot species record vs average distance**

```{r plot sp vs dist,  out.width="65%"}
## get taxa names
key.lab <- unique(mean.dist$class)

## assign point shapes and colors
plot.color <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
point.shape <-  c(0,1,2,3,4,5,6,7,8,9,10)

## plot species richness vs mean distance
p <- ggplot(mean.dist, (aes(x=meanDistance, y=NumbSp, color=class, shape=class))) + theme_bw() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), 
              legend.key.width = unit(0.8, units = "cm"), legend.text = element_text(size=15), 
              legend.key=element_blank()) +
            ylab("Number of species") + 
            xlab("Average distance (m)") +
            labs(title = "Point species richness") +
            geom_point(size = 3, stroke = 0.6) + 
            geom_smooth(method=lm, se=FALSE, fullrange=FALSE, linewidth = 2) +
            scale_color_manual(name = "", labels = key.lab, values = plot.color) +
            scale_shape_manual(name = "", labels = key.lab, values = point.shape)

p

## write plot to a tiff file
tiff("sp_dist.tif", res=300, width = 7, height = 7, unit="in") 
ggplot(mean.dist, (aes(x=meanDistance, y=NumbSp, color=class, shape=class))) + theme_bw() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), 
              legend.key.width = unit(0.8, units = "cm"), legend.text = element_text(size=15), 
              legend.key=element_blank()) +
            ylab("Number of species") + 
            xlab("Average distance (m)") +
            labs(title = "Lumped species richness") +
            geom_point(size = 3, stroke = 0.6) + 
            geom_smooth(method=lm, se=FALSE, fullrange=FALSE, linewidth = 2) +
            scale_color_manual(name = "", labels = key.lab, values = plot.color) +
            scale_shape_manual(name = "", labels = key.lab, values = point.shape)
dev.off()
```

**Get number of conflict sites per province**

```{r conflict freq}
## we are not sure how the frequency of conflicts were scored in the paper as it was not mentioned
## so lets assume that count was done per province 
con.sub.thin.df2 <- con.thin.df
con.sub.thin.df2 <- con.sub.thin.df2 %>% 
                      filter(PROVINCE!="") # filter empty rows

con.freq <- con.sub.thin.df2 %>%
              group_by(PROVINCE) %>%
              count()
## write to a csv
write.csv(con.freq, "conflict.frequency.rows.csv")

## plot number of conflict
ggplot(data=con.freq, aes(x=reorder(PROVINCE, -n), y=n)) +
  geom_bar(stat="identity", fill="#5794db") +
  xlab("Province") +
  ylab("Number of conflict") +
  guides(x =  guide_axis(angle = 75)) 

## combine biodiversity and conflict data (per point)
perpoint.comb <- mean.dist %>% mutate(conflictFreq = if_else(county == "Agusan del Norte", 24,
    if_else(county == "Agusan del Sur", 40,
    if_else(county == "Basilan", 31,
    if_else(county == "Bukidnon", 44,
    if_else(county == "Cagayan de Oro", 0,
    if_else(county == "Camiguin Island", 0,
    if_else(county == "Cotabato", 64,
    if_else(county == "Davao", 0,
    if_else(county == "Davao Oriental", 25,
    if_else(county == "Davao de Oro", 43,
    if_else(county == "Davao del Norte", 19, 
    if_else(county == "Davao del sur", 38,
    if_else(county == "Dinagat Island", 0,
    if_else(county == "General Santos", 0, 
    if_else(county == "Lanao del Norte", 27,
    if_else(county == "Lanao del Sur", 28,
    if_else(county == "Maguindanao", 50, 
    if_else(county == "Misamis Occidental", 10,
    if_else(county == "Misamis Oriental", 16,
    if_else(county == "North Cotabato", 0,
    if_else(county == "Sarangani", 15,
    if_else(county == "South Cotabato", 20,
    if_else(county == "Sultan Kudarat", 29,
    if_else(county == "Sulu", 33,
    if_else(county == "Surigao del Norte", 15,
    if_else(county == "Surigao del Sur", 28,
    if_else(county == "Tawi-Tawi", 12,
    if_else(county == "Zamboanga City", 0,
    if_else(county == "Zamboanga Sibugay", 15,
    if_else(county == "Zamboanga del Norte", 16,
    if_else(county == "Zamboanga del Sur", 34, 0))))))))))))))))))))))))))))))))

## write to a csv
write.csv(perpoint.comb, "perpoint.comb.csv")
```

**Perform Spatial Mixed Models**
```{r spamm}
## perform spatial mixed model with tree cover, tree density and forest canopy height
## these data can be obtained from the links below
## forest canopy height: https://glad.umd.edu/dataset/gedi
## tree density: https://elischolar.library.yale.edu/yale_fes_data/1/
## tree cover: https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html

## read the file with calculated distance from QGIS
bio.con.dist <- read.csv("data/bio.con.dist.csv", header = TRUE, sep = ",")

lon.lat <- subset(bio.con.dist, select = c("decimalLongitude", "decimalLatitude"))

## read raster files
raster.files <- list.files("data/cdh/", pattern = "\\.tif$", full.names = TRUE)

cov.den.hei <- stack(raster.files)

## extract values of tree cover, tree density and forest canopy height

values <- extract(cov.den.hei , lon.lat)

## bind data and get only the required data
bio.con.dist.cdh <- cbind(bio.con.dist, values)

## add conflict frequency data per province
## count conflict frequency per province
con.sub.2 <- con.thin.df
con.sub.2 <- con.sub.2 %>% filter(PROVINCE!="") #filter empty rows

con.freq <- con.sub.2 %>%
              group_by(PROVINCE) %>%
              count()
con.freq

## enter the frequency value manually
bio.con.dist.cdh.cfreq <- bio.con.dist.cdh %>% 
  mutate(conflictFreq = if_else(county == "Agusan del Norte", 24,
    if_else(county == "Agusan del Sur", 40,
    if_else(county == "Basilan", 31,
    if_else(county == "Bukidnon", 44,
    if_else(county == "Cagayan de Oro", 0,
    if_else(county == "Camiguin Island", 0,
    if_else(county == "Cotabato", 64,
    if_else(county == "Davao", 0,
    if_else(county == "Davao Oriental", 25,
    if_else(county == "Davao de Oro", 43,
    if_else(county == "Davao del Norte", 19, 
    if_else(county == "Davao del sur", 38,
    if_else(county == "Dinagat Island", 0,
    if_else(county == "General Santos", 0, 
    if_else(county == "Lanao del Norte", 27,
    if_else(county == "Lanao del Sur", 28,
    if_else(county == "Maguindanao", 50, 
    if_else(county == "Misamis Occidental", 10,
    if_else(county == "Misamis Oriental", 16,
    if_else(county == "North Cotabato", 0,
    if_else(county == "Sarangani", 15,
    if_else(county == "South Cotabato", 20,
    if_else(county == "Sultan Kudarat", 29,
    if_else(county == "Sulu", 33,
    if_else(county == "Surigao del Norte", 15,
    if_else(county == "Surigao del Sur", 28,
    if_else(county == "Tawi-Tawi", 12,
    if_else(county == "Zamboanga City", 0,
    if_else(county == "Zamboanga Sibugay", 15,
    if_else(county == "Zamboanga del Norte", 16,
    if_else(county == "Zamboanga del Sur", 34, 0))))))))))))))))))))))))))))))))

smm.data <- subset(bio.con.dist.cdh.cfreq, select = c(3:8,12, 19,24:27))

## rename year column
colnames(smm.data)[7] <- "Year"

## prepare data for spatial mixed model
## log transform predictors
smm.data$distance.log <- log(smm.data$distance + 1)
smm.data$frequency.log <- log(smm.data$conflictFreq + 1)
smm.data$forestHeight.log <- log(smm.data$forestHeight + 1)
smm.data$treeCov.log <- log(smm.data$treeCov + 1)
smm.data$treeDensity.log <- log(smm.data$treeDensity + 1)

## log transform number of species
smm.data$NumbSp.log <- log(smm.data$NumbSp + 1)

## convert to spatial data
smm.data.sf <- st_as_sf(smm.data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

## remove duplicates
dup.points <- duplicated(st_coordinates(smm.data.sf))
smm.data.sf <- smm.data.sf[!dup.points,]

## remove rows with NAs
smm.data.sf <- na.omit(smm.data.sf)

## create a spatial weights matrix using nearest neighbors
nb.w <- knn2nb(knearneigh(st_coordinates(smm.data.sf), k = 10))
l.w <- nb2listw(nb.w, style = "W") 

plot(nb.w, st_coordinates(smm.data.sf), lwd=.2, col="blue", cex = .5)
title(main = "Spatial Weight Matrix")

## set seed and turn off scientific notation
set.seed(123)
options(scipen = 10)

## assign class and year columns as factors
smm.data.sf$class <- as.factor(smm.data.sf$class)
smm.data.sf$Year <- as.factor(smm.data.sf$Year)

## provide the model equation
eq.mod <- NumbSp.log ~ class + Year + distance.log + frequency.log + forestHeight.log + 
            treeCov.log + treeDensity.log

## fit model
mod <- lm(eq.mod, data=smm.data.sf)

## test regression residuals to examine the residuals of the model using the spatial relationship matrix
## use Morn's correlation
## the null hypothesis states that there is no spatial correlation in the residuals
lm.morantest(mod, l.w)

## or use LaGrange Multiplier
lm.RStests(mod, l.w, test="all")

## stop and use the model results above if there is no spatial correlation (p-value >0.05)
## otherwise, proceed with spatial lag model, spatial error models and etc

## find the best model if necessary
mod <- lm(eq.mod, data = smm.data.sf)

## fit models with fewer predictors
mod1 <- update(mod, . ~ . - class)
mod2 <- update(mod, . ~ . - Year)
mod3 <- update(mod, . ~ . - distance.log)
mod4 <- update(mod, . ~ . - frequency.log)
mod5 <- update(mod, . ~ . - forestHeight.log)
mod6 <- update(mod, . ~ . - treeCov.log)
mod7 <- update(mod, . ~ . - treeDensity.log)

## compare AIC values
aic.values <- AIC(mod, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
aic.values

## identify the model with the lowest AIC
best.model <- rownames(aic.values)[which.min(aic.values$AIC)]
best.model

## use model 2 results to check for spatial correlation
lm.morantest(mod2, l.w)

## or use LaGrange Multiplier
lm.RStests(mod2, l.w, test="all")

## show model 2 summary
## this summary shows some significant results for several classes based on p-value
## however, t-values are too low; thus, these predictors do not contribute much to the model
summary(get(best.model))

## plot results
## get fixed effects coefficients and their confidence intervals
fix.eff <- summary(get(best.model))$coefficients

# create a data frame
fix.eff_df <- as.data.frame(fix.eff)
fix.eff_df$term <- rownames(fix.eff_df)

## plot coefficients
ggplot(fix.eff_df, aes(x = term, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`, 
    ymax = Estimate + 1.96 * `Std. Error`), width = 0.2) +
  theme_minimal() +
  labs(title = "Fixed Effects Estimates", x = "Predictor", y = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## unncomment to reorder plot from lowest to heighest coefficients
##fix.eff_df$term <- factor(fix.eff_df$term, 
##                                levels = fix.eff_df$term[order(fix.eff_df$Estimate)])

##ggplot(fix.eff_df, aes(x = term, y = Estimate)) +
##  geom_point() +
##  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`, 
##    ymax = Estimate + 1.96 * `Std. Error`), width = 0.2) +
##  theme_minimal() +
##  labs(title = "Fixed Effects Estimates", x = "Predictor", y = "Estimate") +
##  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## perform spatial mixed models with land cover data
## land cover data of 1988 was from Swedish Space Corporation (SSC 1988)
## retrieved from 
## https://data.humdata.org/dataset/29a3760f-3170-4555-b5d7-1fbd6cfb5a69?force_layout=desktop
## land cover of 2020 was from the Philippine National Mapping and Resource Information Authority
## data was retrieved from geoportal PH (https://www.geoportal.gov.ph/)
## the national land cover data is ground-thruthed
## this is only for simplistic illustration of sampling artifact
## the original land cover format is shapefile
## we converted it to raster file for fast data processing
## data conversion was done in QGIS using the plugin "rasterize" 
## raster pixel size is 0.0008 or approximately 100 square meters
## all raster files have uniform coordinate reference system, resolution and extent

## read the file with calculated distance from QGIS
bio.con.dist <- read.csv("data/bio.con.dist.csv", header = TRUE, sep = ",")

lon.lat <- subset(bio.con.dist, select = c("decimalLongitude", "decimalLatitude"))

## extract land cover value
values <- extract(lc.data, lon.lat)

## bind data and get only the required data
bio.con.dist.lc <- cbind(bio.con.dist, values)
colnames(bio.con.dist.lc)[24] <- "LC_1988"
colnames(bio.con.dist.lc)[25] <- "LC_2020"

bio.con.dist.lc$LC_1988 <- as.integer(bio.con.dist.lc$LC_1988)
bio.con.dist.lc$LC_2020 <- as.integer(bio.con.dist.lc$LC_2020)

## add land cover descriptions
bio.con.dist.lc.mod <- bio.con.dist.lc %>% 
  mutate(class_1988 = if_else(LC_1988 == 1, "Closed Canopy", 
    if_else(LC_1988 == 2,"Cultivated Area mixed with brushland/grassland",  
    if_else(LC_1988 == 3, "Crop land mixed with coconut plantation",
    if_else(LC_1988 == 4, "Unclassified",
    if_else(LC_1988 == 5, "Open Canopy",
    if_else(LC_1988 == 6, "Arable land, crops mainly cereals and sugar",
    if_else(LC_1988 == 7, "Fishponds derived from mangrove",
    if_else(LC_1988 == 8, "Coral Reef",
    if_else(LC_1988 == 9, "Grassland, grass covering > 70 percent",
    if_else(LC_1988 == 10, "Lake",
    if_else(LC_1988 == 11, "Mangrove vegetation", 
    if_else(LC_1988 == 12, "Built-up Area",
    if_else(LC_1988 == 13, "Quarry",
    if_else(LC_1988 == 14, "Coconut plantations", 
    if_else(LC_1988 == 15, "Siltation pattern in lake",
    if_else(LC_1988 == 16, "Marshy area and swamp",
    if_else(LC_1988 == 17, "Riverbeds",
    if_else(LC_1988 == 18, "Crop land mixed with other plantation",
    if_else(LC_1988 == 19, "Other plantations", "No data"))))))))))))))))))))

bio.con.dist.lc.mod <- bio.con.dist.lc.mod %>% 
  mutate(class_2020 = if_else(LC_2020 == 1, "Grassland", 
    if_else(LC_2020 == 2,"Annual Crop",  
    if_else(LC_2020 == 3, "Marshland/Swamp",
    if_else(LC_2020 == 4, "Perennial Crop",
    if_else(LC_2020 == 5, "Built-up",
    if_else(LC_2020 == 6, "Fishpond",
    if_else(LC_2020 == 7, "Inland Water",
    if_else(LC_2020 == 8, "Mangrove Forest",
    if_else(LC_2020 == 9, "Brush/Shrubs",
    if_else(LC_2020 == 10, "Closed Forest",
    if_else(LC_2020 == 11, "Open Forest", 
    if_else(LC_2020 == 12, "Open/Barren", "No data")))))))))))))

## add conflict frequency data per province
## count conflict frequency per province
con.sub.2 <- con.thin.df
con.sub.2 <- con.sub.2 %>% filter(PROVINCE!="") #filter empty rows

con.freq <- con.sub.2 %>%
              group_by(PROVINCE) %>%
              count()
con.freq

## enter the frequency value manually
bio.con.dist.lc.cfreq <- bio.con.dist.lc.mod %>% 
  mutate(conflictFreq = if_else(county == "Agusan del Norte", 24,
    if_else(county == "Agusan del Sur", 40,
    if_else(county == "Basilan", 31,
    if_else(county == "Bukidnon", 44,
    if_else(county == "Cagayan de Oro", 0,
    if_else(county == "Camiguin Island", 0,
    if_else(county == "Cotabato", 64,
    if_else(county == "Davao", 0,
    if_else(county == "Davao Oriental", 25,
    if_else(county == "Davao de Oro", 43,
    if_else(county == "Davao del Norte", 19, 
    if_else(county == "Davao del sur", 38,
    if_else(county == "Dinagat Island", 0,
    if_else(county == "General Santos", 0, 
    if_else(county == "Lanao del Norte", 27,
    if_else(county == "Lanao del Sur", 28,
    if_else(county == "Maguindanao", 50, 
    if_else(county == "Misamis Occidental", 10,
    if_else(county == "Misamis Oriental", 16,
    if_else(county == "North Cotabato", 0,
    if_else(county == "Sarangani", 15,
    if_else(county == "South Cotabato", 20,
    if_else(county == "Sultan Kudarat", 29,
    if_else(county == "Sulu", 33,
    if_else(county == "Surigao del Norte", 15,
    if_else(county == "Surigao del Sur", 28,
    if_else(county == "Tawi-Tawi", 12,
    if_else(county == "Zamboanga City", 0,
    if_else(county == "Zamboanga Sibugay", 15,
    if_else(county == "Zamboanga del Norte", 16,
    if_else(county == "Zamboanga del Sur", 34, 0))))))))))))))))))))))))))))))))

## subset
smm.data <- subset(bio.con.dist.lc.cfreq, select = c(3:8,12,19,24:28))
colnames(smm.data)[7] <- "Year"

## build a spatial mixed model
## assign land cover, class and year as factors
smm.data$class_1988 <- as.factor(smm.data$class_1988)
smm.data$class_2020 <- as.factor(smm.data$class_2020)
smm.data$class <- as.factor(smm.data$class)
smm.data$Year <- as.factor(smm.data$Year)

## scale distance and frequency data
smm.data$distance.log <- log(smm.data$distance + 1)
smm.data$frequency.log <- log(smm.data$conflictFreq + 1)

## log transform number of species
smm.data$NumbSp.log <- log(smm.data$NumbSp + 1)

## set reference level for land cover
#smm.data$class_2020 <- relevel(smm.data$class_2020, ref = "Built-up")
#smm.data$class_1988 <- relevel(smm.data$class_1988, ref = "Built-up Area")

smm.data$class_2020 <- relevel(smm.data$class_2020, ref = "Closed Forest")
smm.data$class_1988 <- relevel(smm.data$class_1988, ref = "Closed Canopy")

## convert to spatial data
smm.data.sf <- st_as_sf(smm.data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

## remove duplicates
dup.points <- duplicated(st_coordinates(smm.data.sf))
smm.data.sf <- smm.data.sf[!dup.points,]

## create a spatial weights matrix using nearest neighbors
nb.w <- knn2nb(knearneigh(st_coordinates(smm.data.sf), k = 5))
l.w <- nb2listw(nb.w, style = "W") 

#plot(nb.w, st_coordinates(smm.data.sf), lwd=.2, col="blue", cex = .5)
#title(main = "Spatial Weight Matrix")

## set seed and turn off scientific notation
set.seed(123)
options(scipen = 10)

## equation for 1988 land cover
eq.88 <- NumbSp.log ~ class_1988 + class + Year + distance.log + frequency.log 

## fit model for 1988 land cover
mod <- lm(eq.88, data=smm.data.sf)

## test regression residuals to examine the residuals of the model
## use the spatial relationship matrix
lm.morantest(mod, l.w)

## or use LaGrange Multiplier
lm.RStests(mod, l.w, test="all")

## find the best model if necessary
mod <- lm(eq.88, data = smm.data.sf)

## fit models with fewer predictors
mod1 <- update(mod, . ~ . - class_1988)
mod2 <- update(mod, . ~ . - class)
mod3 <- update(mod, . ~ . - Year)
mod4 <- update(mod, . ~ . - distance.log)
mod5 <- update(mod, . ~ . - frequency.log)

## compare AIC values
aic.values <- AIC(mod, mod1, mod2, mod3, mod4, mod5)
aic.values

## identify the model with the lowest AIC
best.model <- rownames(aic.values)[which.min(aic.values$AIC)]
best.model

## use model 3
summary(get(best.model))

## plot results
## get fixed effects coefficients and their confidence intervals
fix.eff <- summary(get(best.model))$coefficients

# create a data frame
fix.eff_df <- as.data.frame(fix.eff)
fix.eff_df$term <- rownames(fix.eff_df)

## plot results
ggplot(fix.eff_df, aes(x = term, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`, 
    ymax = Estimate + 1.96 * `Std. Error`), width = 0.2) +
  theme_minimal() +
  labs(title = "Fixed Effects Estimates with 1988 LC", x = "Predictor", y = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## equation with 2020 land cover
eq.20 <- NumbSp.log ~ class_2020 + class + Year + distance.log + frequency.log

## create a spatial weights matrix using nearest neighbors
nb.w <- knn2nb(knearneigh(st_coordinates(smm.data.sf), k = 10))
l.w <- nb2listw(nb.w, style = "W") 

#plot(nb.w, st_coordinates(smm.data.sf), lwd=.2, col="blue", cex = .5)
#title(main = "Spatial Weight Matrix")

## fit model for 2020 land cover
mod <- lm(eq.20, data=smm.data.sf)

## test regression residuals to examine the residuals of the model
## use the spatial relationship matrix
lm.morantest(mod, l.w)

## or use LaGrange Multiplier
lm.RStests(mod, l.w, test="all")

## find the best model if necessary
mod <- lm(eq.20, data = smm.data.sf)

## fit models with fewer predictors
mod1 <- update(mod, . ~ . - class_2020)
mod2 <- update(mod, . ~ . - class)
mod3 <- update(mod, . ~ . - Year)
mod4 <- update(mod, . ~ . - distance.log)
mod5 <- update(mod, . ~ . - frequency.log)

## compare AIC values
aic.values <- AIC(mod, mod1, mod2, mod3, mod4, mod5)
aic.values

## identify the model with the lowest AIC
best.model <- rownames(aic.values)[which.min(aic.values$AIC)]
best.model

## use model 3
summary(get(best.model))

## plot results
## get fixed effects coefficients and their confidence intervals
fix.eff <- summary(get(best.model))$coefficients

# create a data frame
fix.eff_df <- as.data.frame(fix.eff)
fix.eff_df$term <- rownames(fix.eff_df)

## plot results
ggplot(fix.eff_df, aes(x = term, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = Estimate - 1.96 * `Std. Error`, 
    ymax = Estimate + 1.96 * `Std. Error`), width = 0.2) +
  theme_minimal() +
  labs(title = "Fixed Effects Estimates with 2020 LC", x = "Predictor", y = "Estimate") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

**Calculate proportion of species occurrence records and conflict sites positioned in different land cover types**
```{r land cover, out.width="70%"}
## here we count the number of conflict sites found in different land cover types 
## this demonstrates the artifact of sampling, which was not considered in the paper

## here is the 2020 land cover of Mindanao
knitr::include_graphics("images/lc_2020.png")

## plot data
plot(lc.data$X1988_lc_rast1km)
points(con.thin.df$longitude, con.thin.df$latitude)

con.coor <- cbind(con.thin.df$longitude, con.thin.df$latitude)

## extract data from land cover using conflict data coordinates
con.lc.data <- extract(lc.data, con.coor)
colnames(con.lc.data)[1] <- "LC_1988"
colnames(con.lc.data)[2] <- "LC_2020"

## remove NAs
con.lc.data <- as.data.frame(na.omit(con.lc.data))
con.lc.data$LC_1988 <- as.integer(con.lc.data$LC_1988)
con.lc.data$LC_2020 <- as.integer(con.lc.data$LC_2020)

head(con.lc.data)

## read land cover types code
lc.88.d <- read.csv("data/1988_lc_classmap.csv")
lc.88.d
lc.20.d <- read.csv("data/2020_lc_classmap.csv")
lc.20.d

## add land cover descriptions
con.lc.data.m1 <- con.lc.data %>% 
  mutate(class_1988 = if_else(LC_1988 == 1, "Closed Canopy", 
    if_else(LC_1988 == 2,"Cultivated Area mixed with brushland/grassland",  
    if_else(LC_1988 == 3, "Crop land mixed with coconut plantation",
    if_else(LC_1988 == 4, "Unclassified",
    if_else(LC_1988 == 5, "Open Canopy",
    if_else(LC_1988 == 6, "Arable land, crops mainly cereals and sugar",
    if_else(LC_1988 == 7, "Fishponds derived from mangrove",
    if_else(LC_1988 == 8, "Coral Reef",
    if_else(LC_1988 == 9, "Grassland, grass covering > 70 percent",
    if_else(LC_1988 == 10, "Lake",
    if_else(LC_1988 == 11, "Mangrove vegetation", 
    if_else(LC_1988 == 12, "Built-up Area",
    if_else(LC_1988 == 13, "Quarry",
    if_else(LC_1988 == 14, "Coconut plantations", 
    if_else(LC_1988 == 15, "Siltation pattern in lake",
    if_else(LC_1988 == 16, "Marshy area and swamp",
    if_else(LC_1988 == 17, "Riverbeds",
    if_else(LC_1988 == 18, "Crop land mixed with other plantation",
    if_else(LC_1988 == 19, "Other plantations", "No data"))))))))))))))))))))

con.lc.data.m2 <- con.lc.data %>% 
  mutate(con.lc.data.m1, class_2020 = if_else(LC_2020 == 1, "Grassland", 
    if_else(LC_2020 == 2,"Annual Crop",  
    if_else(LC_2020 == 3, "Marshland/Swamp",
    if_else(LC_2020 == 4, "Perennial Crop",
    if_else(LC_2020 == 5, "Built-up",
    if_else(LC_2020 == 6, "Fishpond",
    if_else(LC_2020 == 7, "Inland Water",
    if_else(LC_2020 == 8, "Mangrove Forest",
    if_else(LC_2020 == 9, "Brush/Shrubs",
    if_else(LC_2020 == 10, "Closed Forest",
    if_else(LC_2020 == 11, "Open Forest", 
    if_else(LC_2020 == 12, "Open/Barren", "No data")))))))))))))

## count the number of conflict sites in each land cover type
## 1988 data
con.lc.count.88 <- con.lc.data.m2 %>%
                      group_by(class_1988) %>%
                      summarise(n())
colnames(con.lc.count.88)[2] <- "count"

con.lc.count.88

ggplot(data=con.lc.count.88, aes(x=reorder(class_1988, -count), y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  labs(title = "Conflict count per LC type (1988)", x = "LC Types", y = "Count") +
  guides(x =  guide_axis(angle = 75)) 

## write to a csv
con.lc.count.88 <- as.data.frame(con.lc.count.88)
write.csv(con.lc.count.88, "con.lc.count.88.csv")

#2020 data
con.lc.count.20 <- con.lc.data.m2 %>%
                      group_by(class_2020) %>%
                      summarise(n())
colnames(con.lc.count.20)[2] <- "count"

con.lc.count.20

ggplot(data=con.lc.count.20, aes(x=reorder(class_2020, -count), y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  labs(title = "Conflict count per LC type (2020)", x = "LC Types", y = "Count") +
  guides(x =  guide_axis(angle = 75)) 

## write to a csv
con.lc.count.20 <- as.data.frame(con.lc.count.20)
write.csv(con.lc.count.20, "con.lc.count.20.csv")

## extract data from land cover using biodiversity data coordinates
bio.coor <- cbind(bio.thin.df$decimalLongitude, bio.thin.df$decimalLatitude)

bio.lc.data <- extract(lc.data, bio.coor)

colnames(bio.lc.data )[1] <- "LC_1988"
colnames(bio.lc.data )[2] <- "LC_2020"

bio.lc.data <- as.data.frame(na.omit(bio.lc.data))

bio.lc.data$LC_1988 <- as.integer(bio.lc.data$LC_1988)
bio.lc.data$LC_2020 <- as.integer(bio.lc.data$LC_2020)

head(bio.lc.data)

## add land cover descriptions
bio.lc.data.m1 <- bio.lc.data %>% 
  mutate(class_1988 = if_else(LC_1988 == 1, "Closed Canopy", 
    if_else(LC_1988 == 2,"Cultivated Area mixed with brushland/grassland",  
    if_else(LC_1988 == 3, "Crop land mixed with coconut plantation",
    if_else(LC_1988 == 4, "Unclassified",
    if_else(LC_1988 == 5, "Open Canopy",
    if_else(LC_1988 == 6, "Arable land, crops mainly cereals and sugar",
    if_else(LC_1988 == 7, "Fishponds derived from mangrove",
    if_else(LC_1988 == 8, "Coral Reef",
    if_else(LC_1988 == 9, "Grassland, grass covering > 70 percent",
    if_else(LC_1988 == 10, "Lake",
    if_else(LC_1988 == 11, "Mangrove vegetation", 
    if_else(LC_1988 == 12, "Built-up Area",
    if_else(LC_1988 == 13, "Quarry",
    if_else(LC_1988 == 14, "Coconut plantations", 
    if_else(LC_1988 == 15, "Siltation pattern in lake",
    if_else(LC_1988 == 16, "Marshy area and swamp",
    if_else(LC_1988 == 17, "Riverbeds",
    if_else(LC_1988 == 18, "Crop land mixed with other plantation",
    if_else(LC_1988 == 19, "Other plantations", "No data"))))))))))))))))))))

bio.lc.data.m2<- bio.lc.data %>% 
  mutate(bio.lc.data.m1, class_2020 = if_else(LC_2020 == 1, "Grassland", 
    if_else(LC_2020 == 2,"Annual Crop",  
    if_else(LC_2020 == 3, "Marshland/Swamp",
    if_else(LC_2020 == 4, "Perennial Crop",
    if_else(LC_2020 == 5, "Built-up",
    if_else(LC_2020 == 6, "Fishpond",
    if_else(LC_2020 == 7, "Inland Water",
    if_else(LC_2020 == 8, "Mangrove Forest",
    if_else(LC_2020 == 9, "Brush/Shrubs",
    if_else(LC_2020 == 10, "Closed Forest",
    if_else(LC_2020 == 11, "Open Forest", 
    if_else(LC_2020 == 12, "Open/Barren", "No data")))))))))))))

## count the number of biodiversity data points in each land cover type
#1988 data
bio.lc.count.88 <- bio.lc.data.m2 %>%
                      group_by(class_1988) %>%
                      summarise(n())
colnames(bio.lc.count.88)[2] <- "count"

bio.lc.count.88

ggplot(data=bio.lc.count.88, aes(x=reorder(class_1988, -count), y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  labs(title = "Biodiversity data count per LC type (1988)", x = "LC Types", y = "Count") +
  guides(x =  guide_axis(angle = 75)) 

bio.lc.count.88 <- as.data.frame(bio.lc.count.88)
write.csv(bio.lc.count.88, "bio.lc.count.88.csv")

#2020 data
bio.lc.count.20 <- bio.lc.data.m2 %>%
                      group_by(class_2020) %>%
                      summarise(n())
colnames(bio.lc.count.20)[2] <- "count"

bio.lc.count.20

bio.lc.count.20 <- as.data.frame(bio.lc.count.20)

ggplot(data=bio.lc.count.20, aes(x=reorder(class_2020, -count), y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  labs(title = "Biodiversity data count per LC type (2020)", x = "LC Types", y = "Count") +
  guides(x =  guide_axis(angle = 75)) 

bio.lc.count.20 <- as.data.frame(bio.lc.count.20)
write.csv(bio.lc.count.20, "bio.lc.count.20.csv")
```