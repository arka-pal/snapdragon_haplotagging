#### Planoles Field Code
#### updated Jun 4, 2021

## Add libraries in a vector 
libraries <- c("tidyverse", "dplyr", "ggplot2")
lapply(libraries, require, character.only=TRUE)

## set working directory 
setwd("/Users/arka/Desktop/Planoles_FWF_data/")

## read files
spatial_dat <- read.csv("20210601/20210601_TrimbleData_cleaned.csv", header=TRUE) #Trimble Data
veg_dat <- read.csv("20210528/20210528_VegTraitData_SS.csv", header=TRUE) #Vegetative Trait Data

## Correcting latitudes and longitudes
#change the NAs to a string 'Center'
spatial_dat$Dir[is.na(spatial_dat$Dir)] <- "Center" 
table(spatial_dat$Dir) # types of directions
#transformation of degrees to radian and proceed to do the cosinus
cos_deg <- function (degrees) 
  {cos(degrees*pi/180)} 
deg2m <- 111319.488 # meters of one latitude degree
spatial_dat$CorrectedLatitude <- ifelse(spatial_dat$Dir == "North", spatial_dat$LatMean + (spatial_dat$HorzCorr/abs(deg2m)), spatial_dat$LatMean) 
spatial_dat$CorrectedLatitude <- ifelse(spatial_dat$Dir == "South", spatial_dat$LatMean - (spatial_dat$HorzCorr/abs(deg2m)), spatial_dat$LatMean)
spatial_dat$CorrectedLongitude <- ifelse(spatial_dat$Dir == "East", spatial_dat$LongMean + (spatial_dat$HorzCorr/abs(deg2m*cos_deg(spatial_dat$LatMean))), spatial_dat$LongMean) # on any other latitude than Equator, 1° of longitude occupies : deg2m*cos(Latitude) metres
spatial_dat$CorrectedLongitude <- ifelse(spatial_dat$Dir == "West", spatial_dat$LongMean - (spatial_dat$HorzCorr/abs(deg2m*cos_deg(spatial_dat$LatMean))), spatial_dat$LongMean)
axis_dist <- function (dis) {dis*sqrt (2)/2} # special formula for inter-cardinal directions
## i haven't done the NW, NE, SW, SE corrections yet. i will upload when i do it. 
#spatial_dat$CorrectedLatitude <- ifelse(spatial_dat$Dir == "NW", spatial_dat$LatMean + (axis_dist(spatial_dat$HorzCorr)/abs(deg2m)), spatial_dat$LatMean)
#spatial_dat$CorrectedLongitude <- ifelse(spatial_dat$Dir == "NW", spatial_dat$LongMean - (axis_dist(spatial_dat$HorzCorr)/abs(deg2m*cos_deg(spatial_dat$LatMean))), spatial_dat$LongMean)


## Merge vegetative and spatial data
#Change plantIDs to lower case for matching
spatial_dat$ID <- tolower(spatial_dat$ID)
spatial_dat$ConfirmID <- tolower(spatial_dat$ConfirmID)
#Merge the two datasets
merged <- full_join(veg_dat, spatial_dat, by='ID')
#write merged file
write_csv(merged, 'merged.csv')


## Plotting for sanity check
ggplot (merged, aes(x=CorrectedLongitude, y=CorrectedLatitude, color=location)) + geom_point() +
  scale_color_manual(values=c("LC" = "orange",
                              "UC" = "orange",
                              "RR" = "orange",
                              "UMF1" = "purple",
                              "MLF1" = "purple",
                              "UYF1" = "yellow",
                              "YLF1" = "yellow",
                              "YUF1" = "yellow",
                              "NA" = "grey",
                              " " ="grey"))
# plot (spatial_dat$CorrectedLongitude, spatial_dat$CorrectedLatitude)
# points (spatial_dat$LongMean, spatial_dat$LatMean, col="black", pch=1)