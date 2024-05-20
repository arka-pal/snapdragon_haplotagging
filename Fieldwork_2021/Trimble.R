### Trimble Data 2021
### written by Arka Pal
### Last updated 16.11.2021

#### new code
setwd("/Users/arka/Desktop/TrimbleDataHZ2021_iter2/")

### read all
carol <- read.csv("Export/all-carol.csv", na.strings = T, header = F, strip.white = T)
crusher <- read.csv("Export/all-crusher.csv", na.strings = T, header = F, strip.white = T)
picard <- read.csv("Export/all-picard.csv", na.strings = T, header = F, strip.white = T)

## ading the extra files
extra_carol <- read.csv("Export/extra-carol.csv", na.strings = T, strip.white = T, header = F)
extra_crusher <- read.csv("Export/extra-crusher.csv", na.strings = T, strip.white = T, header = F)
extra_picard <- read.csv("Export/extra-picard.csv", na.strings = T, strip.white = T, header = F)
FilesOnDevice_carol <- read.csv("Export/FilesOnDevice-carol.csv", na.strings = T, strip.white = T, header = F)
FilesOnDevice_crusher <- read.csv("Export/FilesOnDevice-crusher.csv", na.strings = T, strip.white = T, header = F)

#manipulating 
carol[,1] <- "Carol"
crusher[,1] <- "Crusher"
picard[,1] <- "Picard"
extra_carol[,1] <- "Carol"
extra_crusher[,1] <- "Crusher"
extra_picard[,1] <- "Picard"
FilesOnDevice_carol[,1] <- "Carol"
FilesOnDevice_crusher[,1] <- "Crusher"

#merge 
coords <- rbind (carol, extra_carol, FilesOnDevice_carol,
                 crusher, extra_crusher, FilesOnDevice_crusher,
                 picard, extra_picard)

# name the rows
names(coords) <- c("Source","Latitude","Longitude","Altitude","Plant.ID","Confirm.Plant.ID",
                   "GPS.horz.corr..m.","GPS.alt.corr..m.","Direction.from.GPS","Plant.status",
                   "Num.flower.stems","Num.vegetative.stems","Estim.Flower.No","Tallest.stem.height",
                   "Leaves.to.1st.flower","First.flowered","Perennial","Cut.stems","Total.flowers",
                   "Total.Fuit","Num.aborted.flowers","Comment","veg.flowering","Accessibility",
                   "date.of.measurement","X","X.1","Date","Time","X.2")

# store somecolumns as factors
str(coords)
coords$Plant.ID <- as.factor(coords$Plant.ID)
coords$Confirm.Plant.ID <- as.factor(coords$Confirm.Plant.ID)
coords$Source <- as.factor(coords$Source)
coords$Plant.status <- as.factor(coords$Plant.status)
coords$Direction.from.GPS <- as.factor(coords$Direction.from.GPS)

#delete rows with no plant ID
coords <- coords[c(which(coords$Plant.ID != "" & coords$Confirm.Plant.ID != "")),]
which(coords$Plant.ID == "" & coords$Confirm.Plant.ID =="") #check

#manipulate the dates to order them to know the approx order of sampling
levels(sort(as.factor(coords$Date)))
#change the issues
coords[which(coords$Date == "5/24/2021"),"Date"] <- "05/24/21"
coords[which(coords$Date == "5/26/2021"),"Date"] <- "05/26/21"
coords[which(coords$Date == "5/31/2019"),"Date"] <- "05/31/21"
coords[which(coords$Date == "6/1/2019"),"Date"] <- "06/01/21"
coords[which(coords$Date == "6/16/2021"),"Date"] <- "06/16/21"
coords[which(coords$Date == "6/17/2021"),"Date"] <- "06/17/21"
coords[which(coords$Date == "6/2/2019"),"Date"] <- "06/02/21"
coords[which(coords$Date == "7/19/2021"),"Date"] <- "07/19/21"
coords[which(coords$Date == "7/21/2021"),"Date"] <- "07/21/21"
coords[which(coords$Date == "7/22/2021"),"Date"] <- "07/22/21"
coords[which(coords$Date == "7/23/2021"),"Date"] <- "07/23/21"
coords[which(coords$Date == "7/24/2021"),"Date"] <- "07/24/21"
coords[which(coords$Date == "7/25/2021"),"Date"] <- "07/25/21"
coords[which(coords$Date == "7/26/2021"),"Date"] <- "07/26/21"
coords[which(coords$Date == "7/4/2020"),"Date"] <- "07/04/21"

#manually changed others in excel - errors due to column shift
write.csv(coords,"/Users/arka/Downloads/coords-all.csv", row.names = F)
coords <- read.csv("/Users/arka/Downloads/coords-all.csv",strip.white = T, header = T, na.strings = T)

#check
str(coords)
levels(sort(as.factor(coords$Date)))
coords$Date <- as.factor(coords$Date)

#reorder and put id.generic 
coords <- coords[order(coords$Date, coords$Source),]
coords$ID.generic <- seq(1:nrow(coords))
View(coords)

#matching IDs
coords$match <- 1
coords$match[toupper(coords$Plant.ID) != toupper(coords$Confirm.Plant.ID)] <- 0
unmatched <- coords[which(coords$match == 0),c("Plant.ID","Confirm.Plant.ID")]
nrow(unmatched) # 119 rows - chaned and edited manually in excel. 
write.csv(coords,"/Users/arka/Downloads/coords-all-toMatch.csv", row.names = F)
#unmatched IDs corrected manually


#read file after manual matching and ID correction
coords_raw <- read.csv("/Users/arka/Downloads/coords-All-raw.csv",  header = T, strip.white = T, na.strings = T)
str(coords_raw)
coords_raw$Plant.ID_cap <- toupper(coords_raw$Plant.ID)
#manual cleaning lats and longs 
replicated_IDs <- coords_raw[duplicated(coords_raw$Plant.ID_cap),c("Plant.ID_cap")]

#read file after manual averaging lats and longs ,etc
coords_avg <- read.csv("/Users/arka/Downloads/coords-all-averaged.csv", na.strings = T, strip.white = T, header = T)
str(coords_avg)

#checking which IDs are missing
#lookup files
lookup <- read.csv("/Users/arka/Desktop/TrimbleData-2021/Antirrhinum_LGC_plates_2021_samples.csv",na.strings = T, header=T, strip.white = T)
str(lookup)
lookup$PlantID <- toupper(lookup$PlantID)
names(lookup)[names(lookup) == "PlantID"] <-"Plant.ID_cap"
lookup$confirm_lookup <- lookup$Plant.ID_cap

coords_toMatch <- coords_avg[,c(1:2,6:9,13)]
ID_matching <- merge(lookup, coords_toMatch, by = "Plant.ID_cap", all.x = TRUE, sort = TRUE)
write.csv(ID_matching, "/Users/arka/Downloads/ID_matching.csv", row.names = F)
