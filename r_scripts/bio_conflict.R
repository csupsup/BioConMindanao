---
title: "Supplementary File 2"
date: "20 April 2024"
output:
  pdf_document: default
  html_document:
    df_print: paged
indent: no
---
```{r setup, include=FALSE, message=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
```
This file documents our reanalysis of the dataset used to examine biodiversity and conflict relationships in the Southern Philippines presented in this paper https://doi.org/10.1038/s44185-024-00044-8

**Load required packages**
```{r load package, message=FALSE}
#make sure to install all packages including required dependencies prior to use
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rgdal)
library(raster)
library(lme4)
library(MASS)
library(piecewiseSEM)
```
**Prepare datasets**
```{r read data}
#biodiversity data was obtained from here https://doi.org/10.15468/rtedgk
#and conflict data was from here https://data.humdata.org/dataset/ucdp-data-for-philippines?
#note that only a subset of this dataset (i.e., those from Mindanao and adjacent islands) was used

#biodiversity data
bio <- read.csv("mobios_tab.csv", header = TRUE, sep = ",")
#subset data
bio.sub <- subset(bio, select = c(12, 23, 20, 15, 16))
bio.sub <- bio.sub %>% filter(county!="") #remove rows with no county information
bio.sub <- bio.sub %>% filter(class!="Malacostraca") #exclude this taxon
bio.sub <- bio.sub %>% filter(class!="Bivalvia") #exclude this taxon
bio.sub <- bio.sub %>% filter(class!="Gastropoda")#exclude this taxon
#check province name
#there are rows where the name of the province is uncertain as indicated by "/"
#"Zamboanga¬†del Norte" was changed to "Zamboanga del Norte" prior to import of data
unique(bio.sub$county)
#remove uncertain province names
bio.sub <- bio.sub %>% filter(!grepl('/', county))

head(bio.sub)

#conflict data
con <- read.csv("conflict_tab_89-21.csv", header = TRUE, sep = ",")
con.sub <- subset(con, select = c(1,3,15,18,30,32,33))

head(con.sub)

write.csv(con.sub, "con.sub.csv")
```
**Count recorded species per taxon in each point and within province**
```{r count species, out.width="70%"}
#the methods of how species counts were done in each site was not provided in the paper
#according to the paper analysis was done at the provincial level
#there could be two different ways how the species counts were done (see figure below)
#this part is extremely important which the authors failed to discuss
#the left panel of the figure shows that species counts were aggregated at the provincial level 
#each point (within the province) receives the same value of species count (or species richness)
#the right panel shows that each point has a unique number of species count
#note that these values are crucial for other downstream analyses 
#particularly when relationships of species counts 
#and distance to the nearest conflict sites are considered
#we generated the values for these two scenarios below 
#and test whether any of them would match the results presented in the paper

knitr::include_graphics("bioconflict.png")

#count unique records (i.e., species richness) per point
bio.data.sum <- bio.sub %>%
              group_by(county, class, decimalLatitude, decimalLongitude) %>%
              mutate(NumbSp = n_distinct(scientificName))

head(bio.data.sum)

#export data
write.csv(bio.data.sum, "bio.data.sum.csv")

#count species records per province (this lumps all the data, 
#thus each point will have a single species richness value)
prov.bio.data.sum <- bio.sub %>%
                          group_by(county, class) %>%
                          mutate(NumbSp = n_distinct(scientificName))
head(prov.bio.data.sum)

#export data
write.csv(prov.bio.data.sum, "prov.data.sum.csv")
```
**Calculate the distance of species record to the nearest conflict site/s using QGIS**
```{r cal distance,  out.width="50%"}
#in QGIS, import the "bio.data.sum" and "con.sub" from R as delimited text files. 
#export these files as GeoPackage or Shapefile to convert the data to meters prior to calculation
#use the CRS UTM Zone 51N
knitr::include_graphics("qgis_geo.png")
knitr::include_graphics("qgis_utm.png")
#import back the newly saved data to QGIS
#get the distance of nearest conflict areas using the "Joint attributes by nearest" plugin
#QGIS > Processing Toolbox > Join attributes by nearest
knitr::include_graphics("qgis_join.png")
#distance can be calculated as well without joining the attributes using the "Distance to nearest hub"
#after calculating the distance, a new layer will be created. 
#the new layer contains the calculated distance in the attribute table
#this also includes the x and y coordinates of the nearest site/point
#then export this new layer as a csv file so we can use data in R for other analysis
#repeat the entire process for provincial data (i.e., prov.bio.data.sum)
#figure below shows the nearest conflict point to species record (black line)
#notice the remainder of conflict points are not included
knitr::include_graphics("hub.png")
```

**Calculate the average distance**
```{r cal mean distance}
#species richness per point in each province
#read the file with calculated distance from QGIS
bio.con.dist <- read.csv("bio.con.dist.csv", header = TRUE, sep = ",")
#get average distance by province for each taxon
mean.dist <- bio.con.dist %>%
                group_by(county, class, NumbSp, decimalLatitude, decimalLongitude) %>%
                summarize(meanDistance = mean(distance, na.rm = TRUE))

head(mean.dist)

write.csv(mean.dist, "mean.distance.csv")

#species richness lumped per province
#read the file with calculated distance in QGIS
prov.bio.con.dist <- read.csv("prov.bio.con.dist.csv", header = TRUE, sep = ",")
#get average distance by province for each taxon
prov.mean.dist <- prov.bio.con.dist %>%
                group_by(county, class, NumbSp, decimalLatitude, decimalLongitude) %>%
                summarize(meanDistance = mean(distance, na.rm = TRUE))
write.csv(prov.mean.dist, "prov.mean.distance.csv")

#get average distance by province with considering the long/lat data (i.e., lumped average)
prov.lmp.dist <- prov.bio.con.dist %>%
                group_by(county, class, NumbSp) %>%
                summarize(lmpMeanDistance = mean(distance, na.rm = TRUE))
write.csv(prov.lmp.dist, "prov.lmp.distance.csv")
```
**Plot species record vs average distance**
```{r plot sp vs dist,  out.width="60%"}
key.lab <- unique(mean.dist$class)
plot.color <- c("#7b3dcc", "#cc3d3d", "#000000", "#d970cb", "#8c994d", "#7f404a","#67E3FC")
point.shape <-  c(0,1,2,3,4,5,6,7,8,9,10)

#plot data with species richness per point
p1 <- ggplot(mean.dist, (aes(x=meanDistance, y=NumbSp, color=class, shape=class))) + theme_bw() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), 
              legend.key.width = unit(0.8, units = "cm"), legend.text = element_text(size=15), 
              legend.key=element_blank()) +
            ylab("number of species/record") + 
            xlab("average distance (m)") +
            geom_point(size = 3, stroke = 0.6) + 
            geom_smooth(method=lm, se=FALSE, fullrange=FALSE, size = 2) +
            scale_color_manual(name = "", labels = key.lab, values = plot.color) +
            scale_shape_manual(name = "", labels = key.lab, values = point.shape)

#plot data with lumped species richness
p2 <- ggplot(prov.mean.dist, (aes(x=meanDistance, y=NumbSp, color=class, shape=class))) + theme_bw() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), 
              legend.key.width = unit(0.8, units = "cm"), legend.text = element_text(size=15), 
              legend.key=element_blank()) +
            ylab("number of species/record") + 
            xlab("average distance (m)") +
            geom_point(size = 3, stroke = 0.6) + 
            geom_smooth(method=lm, se=FALSE, fullrange=FALSE, size = 2) +
            scale_color_manual(name = "", labels = key.lab, values = plot.color) +
            scale_shape_manual(name = "", labels = key.lab, values = point.shape)

ggarrange(p1,p2, nrow = 1, ncol = 2, common.legend = TRUE, align = "hv",  widths = 15, heights = 5)

tiff("sp_dist.tif", res=300, width = 7, height = 7, unit="in") 
ggplot(mean.dist, (aes(x=meanDistance, y=NumbSp, color=class, shape=class))) + theme_bw() +
            theme(axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"), 
              legend.key.width = unit(0.8, units = "cm"), legend.text = element_text(size=15), 
              legend.key=element_blank(), legend.position = c(0.83, 0.85), 
              legend.background=element_blank()) +
            ylab("Number of species record") + 
            xlab("Average distance (m)") +
            geom_point(size = 3, stroke = 0.6) + 
            geom_smooth(method=lm, se=FALSE, fullrange=FALSE, size = 2) +
            scale_color_manual(name = "", labels = key.lab, values = plot.color) +
            scale_shape_manual(name = "", labels = key.lab, values = point.shape)
dev.off()
```

**Get number of conflict sites per province**
```{r conflict freq}
#we are not sure how the frequency of conflicts were scored in the paper as it was not mentioned
#so two ways can be considered here
#first, we can get frequency based on the number of unique conflict points within each province
#i.e., per point (unique latitude/longitude)
con.sub.2 <- con.sub
con.sub.2 <- con.sub.2 %>% filter(PROVINCE!="") #filter empty rows

con.freq <- con.sub.2 %>%
              group_by(PROVINCE) %>%
              mutate(ConflicFreq = n_distinct(latitude, longitude)) #use only the unique latitude
con.freq <- con.freq %>% distinct(PROVINCE, .keep_all = TRUE)
write.csv(con.freq, "conflict.frequency.unique.csv")

#second, count each row as a unique report of conflict regardless of coordinates
con.freq2 <- con.sub.2 %>%
              group_by(PROVINCE) %>%
              count()
write.csv(con.freq2, "conflict.frequency.rows.csv")

ggplot(data=con.freq2, aes(x=PROVINCE, y=n)) +
  geom_bar(stat="identity", fill="#5794db") +
  guides(x =  guide_axis(angle = 75)) 

#combine biodiversity and conflict data (lumped)
prov.comb <- prov.lmp.dist %>% mutate(conflictFreq = if_else(county == "Agusan del Norte", 50,
    if_else(county == "Agusan del Sur", 85,
    if_else(county == "Basilan", 329,
    if_else(county == "Bukidnon", 82,
    if_else(county == "Cagayan de Oro", 0,
    if_else(county == "Camiguin Island", 0,
    if_else(county == "Cotabato", 221,
    if_else(county == "Davao", 0,
    if_else(county == "Davao Oriental", 52,
    if_else(county == "Davao de Oro", 113,
    if_else(county == "Davao del Norte", 56, 
    if_else(county == "Davao del sur", 181,
    if_else(county == "Dinagat Island", 0,
    if_else(county == "General Santos", 0, 
    if_else(county == "Lanao del Norte", 66,
    if_else(county == "Lanao del Sur", 181,
    if_else(county == "Maguindanao", 352, 
    if_else(county == "Misamis Occidental", 15,
    if_else(county == "Misamis Oriental", 38,
    if_else(county == "North Cotabato", 0,
    if_else(county == "Sarangani", 25,
    if_else(county == "South Cotabato", 55,
    if_else(county == "Sultan Kudarat", 62,
    if_else(county == "Sulu", 393,
    if_else(county == "Surigao del Norte", 29,
    if_else(county == "Surigao del Sur", 72,
    if_else(county == "Tawi-Tawi", 17,
    if_else(county == "Zamboanga City", 0,
    if_else(county == "Zamboanga Sibugay", 18,
    if_else(county == "Zamboanga del Norte", 34,
    if_else(county == "Zamboanga del Sur", 106, 0))))))))))))))))))))))))))))))))
write.csv(prov.comb, "prov.lmp.comb.csv")

#combine biodiversity and conflict data (per point)
perpoint.comb <- mean.dist %>% mutate(conflictFreq = if_else(county == "Agusan del Norte", 50,
    if_else(county == "Agusan del Sur", 85,
    if_else(county == "Basilan", 329,
    if_else(county == "Bukidnon", 82,
    if_else(county == "Cagayan de Oro", 0,
    if_else(county == "Camiguin Island", 0,
    if_else(county == "Cotabato", 221,
    if_else(county == "Davao", 0,
    if_else(county == "Davao Oriental", 52,
    if_else(county == "Davao de Oro", 113,
    if_else(county == "Davao del Norte", 56, 
    if_else(county == "Davao del sur", 181,
    if_else(county == "Dinagat Island", 0,
    if_else(county == "General Santos", 0, 
    if_else(county == "Lanao del Norte", 66,
    if_else(county == "Lanao del Sur", 181,
    if_else(county == "Maguindanao", 352, 
    if_else(county == "Misamis Occidental", 15,
    if_else(county == "Misamis Oriental", 38,
    if_else(county == "North Cotabato", 0,
    if_else(county == "Sarangani", 25,
    if_else(county == "South Cotabato", 55,
    if_else(county == "Sultan Kudarat", 62,
    if_else(county == "Sulu", 393,
    if_else(county == "Surigao del Norte", 29,
    if_else(county == "Surigao del Sur", 72,
    if_else(county == "Tawi-Tawi", 17,
    if_else(county == "Zamboanga City", 0,
    if_else(county == "Zamboanga Sibugay", 18,
    if_else(county == "Zamboanga del Norte", 34,
    if_else(county == "Zamboanga del Sur", 106, 0))))))))))))))))))))))))))))))))
write.csv(perpoint.comb, "perpoint.comb.csv")
```

**Perform GLM**
```{r glm}
#glm using lumped species richness per province
#assign taxa as a factor
prov.comb$class <- factor(prov.comb$class)

#plot data
h1 <- ggplot(prov.comb, aes(x = lmpMeanDistance, y = NumbSp)) + 
        geom_point(pch = 1) + geom_smooth(method = lm)
h2 <- ggplot(prov.comb, aes(x = conflictFreq, y = NumbSp)) + 
        geom_point(pch = 1) + geom_smooth(method = lm)
ggarrange(h1, h2, nrow = 1, ncol = 2, align = "hv")

#perform glm
glm.res.1 <- glm(NumbSp ~ conflictFreq + lmpMeanDistance + class, 
                 family="poisson"(link="log"), data = prov.comb)
glm.res.2 <- glm(NumbSp ~ conflictFreq + class, 
                 family="poisson"(link="log"), data = prov.comb)
glm.res.3 <- glm(NumbSp ~ lmpMeanDistance + class, 
                 family="poisson"(link="log"), data = prov.comb)
glm.res.4 <- glm(NumbSp ~ conflictFreq + lmpMeanDistance, 
                 family="poisson"(link="log"), data = prov.comb)

#compare glm results
AIC(glm.res.1, glm.res.2, glm.res.3, glm.res.4)
stepAIC(glm.res.1, direction='both') #use stepAIC function

#get summary for best model
summary(glm.res.1)

#glm with taxa as a random variable
pred.scale <- scale(prov.comb[4:5]) #scale predictors
prov.comb.s <- cbind(prov.comb[2:3], pred.scale)

glm.res.5 <- glmer(NumbSp ~ conflictFreq + lmpMeanDistance + (1|class), 
                 family = "poisson"(link ="log"), data = prov.comb.s)
summary(glm.res.5)

#get R^2
rsquared(glm.res.5, method="trigamma")

#use data per point
#assign taxa as a factor
perpoint.comb$class <- factor(perpoint.comb$class)

#plot data
h1 <- ggplot(perpoint.comb, aes(x = meanDistance, y = NumbSp)) + 
        geom_point(pch = 1) + geom_smooth(method = lm)
h2 <- ggplot(perpoint.comb, aes(x = conflictFreq, y = NumbSp)) + 
        geom_point(pch = 1) + geom_smooth(method = lm)
ggarrange(h1, h2, nrow = 1, ncol = 2, align = "hv")

#perform glm
glm.res.1 <- glm(NumbSp ~ conflictFreq + meanDistance + class, 
                 family = "poisson"(link ="log"), data = perpoint.comb)
glm.res.2 <- glm(NumbSp ~ conflictFreq + class, 
                 family = "poisson"(link ="log"), data = perpoint.comb)
glm.res.3 <- glm(NumbSp ~ meanDistance + class, 
                 family = "poisson"(link ="log"), data = perpoint.comb)
glm.res.4 <- glm(NumbSp ~ conflictFreq + meanDistance, 
                 family = "poisson"(link ="log"), data = perpoint.comb)

#compare glm results
AIC(glm.res.1, glm.res.2, glm.res.3, glm.res.4)
stepAIC(glm.res.1, direction = 'both') #use stepAIC function

#get summary for best model
summary(glm.res.1)

#glm with taxa as a random variable
pred.scale <- scale(perpoint.comb[6:7]) #scale predictors
perpoint.comb.s <- cbind(perpoint.comb[2:3], pred.scale)

glm.res.5 <- glmer(NumbSp ~ conflictFreq + meanDistance + (1|class), 
                 family = "poisson"(link = "log"), data = perpoint.comb.s)
summary(glm.res.5)

#get R^2
rsquared(glm.res.5, method = "trigamma")
```

**Calculate the proportion of species occurrence records and conflict sites positioned in different land cover types**
```{r land cover, out.width="70%"}
#here we calculate the number of sites found in different land cover types 
#this demonstrates the artifact of sampling, which was not considered in the paper
#land cover data of 1988 was from Swedish Space Corporation (SSC 1988)
#land cover of 2020 was from the Philippine National Mapping and Resource Information Authority
#data was retrieved from geoportal PH (https://www.geoportal.gov.ph/)
#we used land cover as a proxy for tree density and forest canopy data 
#the national land cover data is ground-thruthed
#this is only for simplistic illustration of sampling artifact
#the original land cover format is shapefile
#we converted it to raster file for fast data processing
#data conversion was done in QGIS using the plugin "rasterize" 
#raster pixel size is 0.0008 or approximately 100 meters
#all raster files have the uniform coordinate reference system, resolution and extent

#here is the 2020 land cover of Mindanao
knitr::include_graphics("lc_2020.png")

#read data
files <- list.files("land cover/", pattern = "tif", full.names=TRUE)
raster.files <- lapply(files, raster)

raster.files

#resample raster files to have a uniform dimension/resolution
standard <- raster.files[[1]]
raster.resamp <- list(standard)
for (i in 2:length(raster.files)) {
  raster.resamp[[i]] <- resample(raster.files[[i]], standard,
                                 method='bilinear')}
#stack raster files
lc <- stack(raster.files)

#plot data
plot(lc$X1988_lc_mindanao_rast100m)
points(con.sub$longitude, con.sub$latitude)

con.coor <- cbind(con.sub$longitude, con.sub$latitude)

#extract data from land cover using conflict data coordinates
con.lc.data <- extract(lc, con.coor)
colnames(con.lc.data)[1] <- "LC_1988"
colnames(con.lc.data)[2] <- "LC_2020"
con.lc.data <- as.data.frame(na.omit(con.lc.data))
head(con.lc.data)

#land cover types code
lc.88.d <- read.csv("1988_lc_mindanao_code.csv")
lc.88.d
lc.20.d <- read.csv("2020_lc_mindanao_code.csv")
lc.20.d

#add land cover descriptions
con.lc.data.m1 <- con.lc.data %>% mutate(class_1988 = if_else(LC_1988 == 4, "Lake", 
    if_else(LC_1988 == 5,"Cultivated Area mixed with brushland/grassland",  
    if_else(LC_1988 == 7, "Built-up Area",
    if_else(LC_1988 == 13, "Quarry",
    if_else(LC_1988 == 15, "Unclassified",
    if_else(LC_1988 == 37, "Coconut plantations",
    if_else(LC_1988 == 42, "Open Canopy",
    if_else(LC_1988 == 154, "Crop land mixed with coconut plantation",
    if_else(LC_1988 == 487, "Arable land, crops mainly cereals and sugar",
    if_else(LC_1988 == 504, "Coral Reef",
    if_else(LC_1988 == 557, "Riverbeds", 
    if_else(LC_1988 == 667, "Open/Siltation pattern in lake",
    if_else(LC_1988 == 683, "Mangrove vegetation",
    if_else(LC_1988 == 706, "Fishponds derived from mangrove", 
    if_else(LC_1988 == 728, "Marshy area and swamp",
    if_else(LC_1988 == 746, "Closed Canopy",
    if_else(LC_1988 == 756, "Crop land mixed with other plantation",
    if_else(LC_1988 == 776, "Other plantations",
    if_else(LC_1988 == 840, "Grassland, grass covering > 70 percent",
    if_else(LC_1988 == 867, "Open Canopy", "No data")))))))))))))))))))))

con.lc.data.m2 <- con.lc.data %>% mutate(con.lc.data.m1, class_2020 = if_else(LC_2020 == 1, "Grassland", 
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

#count the number of conflict sites in each land cover type
#1988 data
con.lc.count.88 <- con.lc.data.m2 %>%
                      group_by(class_1988) %>%
                      summarise(n())
colnames(con.lc.count.88)[2] <- "count"

ggplot(data=con.lc.count.88, aes(x=class_1988, y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  guides(x =  guide_axis(angle = 75)) 

con.lc.count.88 <- as.data.frame(con.lc.count.88)
write.csv(con.lc.count.88, "con.lc.count.88.csv")

#2020 data
con.lc.count.20 <- con.lc.data.m2 %>%
                      group_by(class_2020) %>%
                      summarise(n())
colnames(con.lc.count.20)[2] <- "count"

ggplot(data=con.lc.count.20, aes(x=class_2020, y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  guides(x =  guide_axis(angle = 75)) 

con.lc.count.20 <- as.data.frame(con.lc.count.20)
write.csv(con.lc.count.20, "con.lc.count.20.csv")

#extract data from land cover using biodiversity data coordinates
bio.coor <- cbind(bio.sub$decimalLongitude, bio.sub$decimalLatitude)
bio.lc.data <- extract(lc, bio.coor)
colnames(bio.lc.data )[1] <- "LC_1988"
colnames(bio.lc.data )[2] <- "LC_2020"
bio.lc.data <- as.data.frame(na.omit(bio.lc.data ))
head(bio.lc.data)

#add land cover descriptions
bio.lc.data.m1 <- bio.lc.data %>% mutate(class_1988 = if_else(LC_1988 == 4, "Lake", 
    if_else(LC_1988 == 5,"Cultivated Area mixed with brushland/grassland",  
    if_else(LC_1988 == 7, "Built-up Area",
    if_else(LC_1988 == 13, "Quarry",
    if_else(LC_1988 == 15, "Unclassified",
    if_else(LC_1988 == 37, "Coconut plantations",
    if_else(LC_1988 == 42, "Open Canopy",
    if_else(LC_1988 == 154, "Crop land mixed with coconut plantation",
    if_else(LC_1988 == 487, "Arable land, crops mainly cereals and sugar",
    if_else(LC_1988 == 504, "Coral Reef",
    if_else(LC_1988 == 557, "Riverbeds", 
    if_else(LC_1988 == 667, "Open/Siltation pattern in lake",
    if_else(LC_1988 == 683, "Mangrove vegetation",
    if_else(LC_1988 == 706, "Fishponds derived from mangrove", 
    if_else(LC_1988 == 728, "Marshy area and swamp",
    if_else(LC_1988 == 746, "Closed Canopy",
    if_else(LC_1988 == 756, "Crop land mixed with other plantation",
    if_else(LC_1988 == 776, "Other plantations",
    if_else(LC_1988 == 840, "Grassland, grass covering > 70 percent",
    if_else(LC_1988 == 867, "Open Canopy", "No data")))))))))))))))))))))

bio.lc.data.m2<- bio.lc.data %>% mutate(bio.lc.data.m1, class_2020 = if_else(LC_2020 == 1, "Grassland", 
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

#count the number of conflict sites in each land cover type
#1988 data
bio.lc.count.88 <- bio.lc.data.m2 %>%
                      group_by(class_1988) %>%
                      summarise(n())
colnames(bio.lc.count.88)[2] <- "count"

ggplot(data=bio.lc.count.88, aes(x=class_1988, y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  guides(x =  guide_axis(angle = 75)) 

bio.lc.count.88 <- as.data.frame(bio.lc.count.88)
write.csv(bio.lc.count.88, "bio.lc.count.88.csv")

#2020 data
bio.lc.count.20 <- bio.lc.data.m2 %>%
                      group_by(class_2020) %>%
                      summarise(n())
colnames(bio.lc.count.20)[2] <- "count"

bio.lc.count.20 <- as.data.frame(bio.lc.count.20)

ggplot(data=bio.lc.count.20, aes(x=class_2020, y=count)) +
  geom_bar(stat="identity", fill="#5794db") +
  guides(x =  guide_axis(angle = 75)) 

bio.lc.count.20 <- as.data.frame(bio.lc.count.20)
write.csv(bio.lc.count.20, "bio.lc.count.20.csv")
```
