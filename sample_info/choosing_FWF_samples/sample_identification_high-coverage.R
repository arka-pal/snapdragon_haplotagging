fwf <- read.csv("haplotag_samples_new.csv", header = T, na.strings = T, strip.white = T)
str(fwf)

fwf$lab <- paste(fwf$ID,'',fwf$Red,'',fwf$Yellow)
fwf$Red <- as.numeric(fwf$Red)
fwf$Yellow <- as.numeric(fwf$Yellow)

fwf$Red[321] <- 100 #NA value coerced with 100
fwf$Yellow[321] <- 100 #NA value coerced with 100

for (i in c(1:nrow(fwf))){
  fwf$c[i] <- if (fwf$Red[i] <= 0.5 && fwf$Yellow[i] >= 2.5){"Yellow"}
  else if(fwf$Red[i] > 3 && fwf$Yellow[i] <= 1){"Magenta"}
  else if(fwf$Red[i] >= 1.5 && fwf$Red[i] <= 2.5 && fwf$Yellow[i] <= 1){"Pink"}
  else if(fwf$Red[i] > 2.5 && fwf$Yellow[i] > 2){"Orange"}
  else if(fwf$Red[i] <= 1 && fwf$Yellow[i] <= 1){"lightgrey"}
  else{"Black"}
  }

p <- ggplot(fwf, aes(label=lab)) + geom_point(colour=fwf$c,aes(AvgLong, AvgLat)) 
ggplotly(p) 
