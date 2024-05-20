library(plotly)
library(ggplot2)

setwd("/Volumes/apal/from_mac_202111123_desktop/TrimbleDataCleaning-2021/sample_identification_for_FWF")

all_coords <- read.csv("coords_ALL-2021-v3.1_20211117.csv", strip.white = T, na.strings = T, header = T)
str(all_coords)
all_coords$AvgLat <- as.numeric(all_coords$AvgLat)
all_coords$AvgLong <- as.numeric(all_coords$AvgLong)
all_Alive <- all_coords[which(all_coords$Plant.status == "Alive"),]
all_Alive <- all_Alive[which(all_Alive$Flag == 0),]

colour_scores <- read.csv("/Users/arka/Downloads/20210801_HZphotoscores2021_AP.csv", header = T)
str(colour_scores)
colour_scores <- colour_scores[,-c(11:19)]
colour_scores$Plant.ID_cap <- toupper(colour_scores$Plant_No)
colour_scores$Red <- as.numeric(colour_scores$Red)
colour_scores$Yellow <- as.numeric(colour_scores$Yellow)

all_scored <- merge(colour_scores, all_Alive, by = "Plant.ID_cap", all.x = T)
str(all_scored)

## picking individuals from different clusters 

#MF_far
MF2 <- all_scored[which(all_scored$AvgLong <= 2.0964 & all_scored$AvgLong >= 2.0953 & all_scored$AvgLat <= 42.3210 & all_scored$AvgLat >= 42.3207),]
nrow(MF2) #59
MF2 <- MF2[order(MF2$AvgLong, decreasing = F),]
tmp_MF2 <- data.frame(ID = c(MF2$Plant.ID_cap), Order=c(1:nrow(MF2)), Location=rep("MF_far",nrow(MF2)))

#MF_near
# MF1 <- all_scored[which(all_scored$AvgLong <= 2.0920 & all_scored$AvgLong >= 2.0903 & all_scored$AvgLat <= 42.3224 & all_scored$AvgLat >= 42.3220),]
# nrow(MF1)#133
MF1 <- all_scored[which(all_scored$AvgLong <= 2.0894 & all_scored$AvgLong >= 2.0885 & all_scored$AvgLat <= 42.3221 & all_scored$AvgLat >= 42.3217),]
nrow(MF1)#71
MF1 <- MF1[order(MF1$AvgLong, decreasing = F),]
tmp_MF1 <- data.frame(ID = c(MF1$Plant.ID_cap), Order=c(1:nrow(MF1)), Location=rep("MF_near",nrow(MF1)))

#YF_near
YF1 <- all_scored[which(all_scored$AvgLong <= 2.0636 & all_scored$AvgLong >= 2.061 & all_scored$AvgLat <= 42.326 & all_scored$AvgLat >= 42.324),]
YF1 <- YF1[order(YF1$AvgLong, decreasing = F),]
nrow(YF1) #40
tmp_YF1 <- data.frame(ID = c(YF1$Plant.ID_cap), Order=c(1:nrow(YF1)), Location=rep("YF_near",nrow(YF1)))

#YF_far
YF2 <- (all_scored[which(all_scored$AvgLong <= 2.0526 & all_scored$AvgLong >= 2.052 & all_scored$AvgLat <= 42.3257 & all_scored$AvgLat >= 42.3252),])
nrow(YF2) #75
YF2 <- YF2[order(YF2$AvgLong, decreasing = F),]
tmp_YF2 <- data.frame(ID = c(YF2$Plant.ID_cap), Order=c(1:nrow(YF2)), Location=rep("YF_far",nrow(YF2)))

#combine all IDs
flank_list <- rbind(tmp_YF2, tmp_YF1, tmp_MF1, tmp_MF2)
write.csv(flank_list,"list.csv",row.names = F)


coords <- ggplot(all_Alive) +
  geom_point(aes(AvgLong, AvgLat)) +
  geom_point(colour="red",aes(AvgLong, AvgLat), data = MF1) +
  geom_point(colour="red",aes(AvgLong, AvgLat), data = MF2) +
  geom_point(colour="yellow",aes(AvgLong, AvgLat), data = YF1) +
  geom_point(colour="yellow",aes(AvgLong, AvgLat), data = YF2) 
  #+ geom_point(colour="green",aes(AvgLong, AvgLat), data = reds_in_YF[which(reds_in_YF$Red == 4 & reds_in_YF$Yellow < 1.5),])

coords
ggplotly(coords)

par(mfrow=c(2,2))
plot(jitter(Yellow)~jitter(Red), data = YF2, main = "YF_far", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")
plot(jitter(Yellow)~jitter(Red), data = YF1, main = "YF_near", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")
plot(jitter(Yellow)~jitter(Red), data = MF1, main = "MF_near", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")
plot(jitter(Yellow)~jitter(Red), data = MF2, main = "MF_far", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")

## picking outlier individuals
reds_in_YF <- all_scored[which(all_scored$AvgLong <= 2.0602),]
yellows_in_MF <- all_scored[which(all_scored$AvgLong >= 2.0903),]

#Zxxxx
old <- read.csv("/Users/arka/Downloads/locatpheno_filtered_AP.csv", header = T)
str(old)
old_MF <- old[which(old$Location > 13700),]
old_YF <- old[which(old$Location < 13000),]

#PAxxxx
coords_2020 <- read.csv("/Users/arka/Downloads/20200304_AllGPS_PipelineReady_CB.csv", header = T)
colours_2020 <- read.csv("/Users/arka/Downloads/HZPhotoScores2020_AP.csv", header = T)
colours_2020$Plant.ID_cap <- toupper(colours_2020$Plant.No)
merged_2020 <- merge(colours_2020, coords_2020, by="Plant.ID_cap", all = T)
YF_2020 <- merged_2020[which(merged_2020$Longitude <= 2.0602),]
MF_2020 <- merged_2020[which(merged_2020$Longitude >= 2.0903),]

#plot old and new data
par(mfrow=c(3,2))
plot(jitter(Yellow)~jitter(Red), data = reds_in_YF, main = "Yellow Flank (Upper + Lower Road)", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")
plot(jitter(Yellow)~jitter(Red), data = yellows_in_MF, main = "Magenta Flank (Upper + Lower Road)", xlim=c(0,5), ylim=c(0,4), xlab="Red Score", ylab="Yellow Score")
plot (jitter(Yellow_final)~jitter(Red_final), old_YF, xlim=c(0,5), ylim=c(0,4), main="Zxxxx Yellow Flank (Upper + Lower Road)")
plot (jitter(Yellow_final)~jitter(Red_final), old_MF, xlim=c(0,5), ylim=c(0,4), main="Zxxxx Magenta Flank (Upper + Lower Road)")
plot (jitter(Yellow)~jitter(Red), YF_2020, xlim=c(0,5), ylim=c(0,4), main="PAxxxx Yellow Flank (Upper + Lower Road)")
plot (jitter(Yellow)~jitter(Red), MF_2020, xlim=c(0,5), ylim=c(0,4), main="PAxxxx Magenta Flank (Upper + Lower Road)")


#reds in YF
reds_in_YF[which(reds_in_YF$Red >= 4 & reds_in_YF$Yellow <= 1.5),]$Plant.ID_cap
reds_in_YF[which(reds_in_YF$Red >= 4 & reds_in_YF$Yellow <= 1.5),]
#yellows in MF
yellows_in_MF[which(yellows_in_MF$Red < 1 & yellows_in_MF$Yellow >= 1.5),]$Plant.ID_cap
yellows_in_MF[which(yellows_in_MF$Red < 1 & yellows_in_MF$Yellow >= 1.5),]$Plant.ID_cap
#reds in YF from Zxxxx
old_YF[which(old_YF$Red_final >= 4 & old_YF$Yellow_final <= 1.5),]$PlantID_final
old_YF[which(old_YF$Red_final >= 4 & old_YF$Yellow_final <= 1.5),]
#yellows in MF from Zxxxx
old_MF[which(old_MF$Red_final < 1 & old_MF$Yellow_final >= 1.5),]$PlantID_final
old_MF[which(old_MF$Red_final < 1 & old_MF$Yellow_final >= 1.5),]
#reds in YF from PAxxxx
YF_2020[which(YF_2020$Red >= 4 & YF_2020$Yellow <= 1.5),]$Plant.ID_cap
YF_2020[which(YF_2020$Red >= 4 & YF_2020$Yellow <= 1.5),]
#yellows in MF from PAxxxx
MF_2020[which(MF_2020$Red < 1 & MF_2020$Yellow >= 1.5),]$Plant.ID_cap
MF_2020[which(MF_2020$Red < 1 & MF_2020$Yellow >= 1.5),]
