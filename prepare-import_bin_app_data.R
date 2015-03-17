#
#      Read, import and bin appendicularians (for Carbon content)
#
#--------------------------------------------------------------------------#

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- TRUE

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

## { Import data ----------------------------------------------------------
bioData <- "raw_biological_data"

AFiles <- list.files(bioData, pattern=glob2rx("*AppsCarbon*csv"), full=TRUE)
apps <- adply(AFiles, 1, read.csv, stringsAsFactors=FALSE)

# identify transects (based on the day number)
apps$X1 <- AFiles[apps$X1]
apps$X1 <- str_replace(apps$X1, bioData, "")
apps$transect <- as.numeric(str_sub(apps$X1, 4, 5)) - 14

# create a true date+time column
apps$day <- as.numeric(str_sub(apps$X1, 4, 5))

# shift by one day when we cross midnight
apps$day <- ifelse(apps$hour >= 18 & apps$hour <= 23, apps$day, apps$day+1)

# convert that into POSIXct
apps$dateTime <- as.POSIXct(str_c("2010-10-", apps$day, " ", sprintf("%02i",apps$hour), ":", sprintf("%02i",apps$min), ":", sprintf("%02i",apps$sec), ".", apps$s.1000), tz="America/Los_Angeles")

apps <- apps[ , ! names(apps) %in% c("X1", "day", "hour", "min", "sec", "s.1000", "in.hse", "out.hse")]
apps <- rename(apps, replace=c("kow"="kowalevskiid", "fri"="fritillarid"))

appsT <- melt(apps, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")

appsT$group <- "Tunicates"
appsT$sub <- 20

appsT$dateTime <- appsT$dateTime - 3 * 3600

# }

##{ Read in physical data and bin (same code as in prepare-bin_by_depth.R)--------------------

phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="America/Los_Angeles")

# depth bin size in m
binSize <- 1

maxD <- max(phy$depth, na.rm=T)

bins=seq(0, maxD, by=binSize)

# per cast, interpolate the times corresponding to the crossing of the depth bins
timeBins <- dlply(phy, ~transect+cast+down.up, function(x, bins) {  
  # compute corresponding crossing times
  timeBins <- na.omit(approx(x$depth, x$dateTime, bins, rule=1)$y)
  return(timeBins)
}, bins=bins, .progress="text")
# combine the elements of the lists into just one vector
timeBins <- unique(sort(do.call(c, timeBins)))

# add times at which we shift between up and downcasts
d.u <- ifelse(phy$down.up == "down", 1, 2)
d.u <- diff(d.u)
changes <- which(d.u != 0)
changes <- data.frame(phy$dateTime[changes], phy$dateTime[changes+1])
changeTimes <- aaply(changes, 1, function(x) {
  mean(unlist(x))
}, .expand=FALSE)

timeBins <- sort(c(timeBins, changeTimes))

# make it a proper time
timeBins <- as.POSIXct(timeBins, origin="1970-01-01 00:00:00.00 GMT", tz="America/Los_Angeles")

# verify where the cuts are
dd <- data.frame(dateTime=timeBins)
# find the depth/long/transect number where the cuts will happen
dd$depth <- approx(phy$dateTime, phy$depth, timeBins)$y
dd$long <- approx(phy$dateTime, phy$long, timeBins)$y
dd$transect <- approx(phy$dateTime, phy$transect, timeBins)$y

## Bin physical Data ##

# create the bins (and ensure labels are simple, consistent and unique)
phy$dateTimeB <- cut(phy$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# convert from factor to numeric
phy$dateTimeB <- as.numeric(phy$dateTimeB)

# calculate the time duration of each bin
# assign the timeBin values to the phy dataframe, convert to POSIX
phy$timeBins <- cut(phy$dateTime, breaks=timeBins, labels=timeBins[1:length(timeBins)-1])
phy$timeBins <- as.POSIXct(phy$timeBins, origin="1970-01-01 00:00:00.00 GMT", tz="America/Los_Angeles")

# cut by transect, because the time bin duration values need to be calculated per bin
binDuration <- dlply(phy, ~transect, function(x){
  # because time bin values are repeated, find just the unique values
  uniqtimeBins <- unique(x$timeBins)
  
  # calculate the duration of the first and last bin
  durationFirstBin <- difftime(uniqtimeBins[2], min(x$dateTime), units="secs")
  durationLastBin <- difftime(max(x$dateTime), max(uniqtimeBins, na.rm=T), units="secs")
  
  # calculate the time difference between bins
  binDuration <- diff(uniqtimeBins)
  
  # add in the duration of the first bin and the last bin
  binDuration <- c(durationFirstBin, binDuration[2:length(binDuration)], durationLastBin)
  
  return(binDuration)
})

# convert list to data frame, then combine the  
binDuration <- rbind(data.frame(a=as.vector(binDuration$`1`)), data.frame(a=as.vector(binDuration$`2`)), data.frame(a=as.vector(binDuration$`3`)))
binDuration <- binDuration[complete.cases(binDuration),]
binDuration <- as.difftime(binDuration, units="secs")

binDuration <- data.frame(duration=binDuration, dateTimeB=as.numeric(unique(phy$dateTimeB)))
binDuration[1,2] <- 0

# remove timeBins column
phy <- phy[, which(names(phy) != "timeBins")]

# call the first datetimeB "0" instead of NA
phy[which(is.na(phy[phy$transect==1 & phy$cast==1,"dateTimeB"])),"dateTimeB"] <- 0

phy[which(is.na(phy$dateTimeB)),"dateTimeB"] <- as.numeric(max(phy$dateTimeB, na.rm=T))

# add another dateTimeB column
phy$dateTimeB2 <- phy$dateTimeB

# compute average properties per bin
phyB <- ddply(phy, ~transect+cast+front+down.up+dateTimeB, function(x) {
  
  # means of variables inside the bin
  out <- as.data.frame(as.list(colMeans(x[,c("depth", "lat", "long", "temp", "salinity", "swRho", "fluoro", "oxygen", "irradiance", "vertical.vel", "pitch", "velocity")], na.rm=T)))
  
  # retrieve computed duration
  index <- as.numeric(unique(x$dateTimeB2))
  
  out$duration <- as.numeric(binDuration[which(binDuration$dateTimeB == index), "duration"])
  
  if (out$duration > 60){
    out$duration <-  as.numeric(diff(range(x$dateTime, na.rm=T)))
  }
  
  out$n <- nrow(x)
  
  # compute average date
  out$dateTime <- mean(x$dateTime, na.rm=T)
  
  # when horizontal velocity is not available, assign something to it to be able to compte actual velocity
  if (is.na(out$velocity)) {
    out$velocity <- 2.5119
  }
  
  # compute volume sampled in that bin (in m^3)
  #             ISIIS field in m^2  *  speed in m/s           *  duration in s  
  out$volume <- 0.13 * 0.50         *  (out$velocity)         *  out$duration
  
  return(out)
}, .progress="text", .parallel=parallel)

# }

## { Bin apps data & Join with physical data ------------------------------------------------

# create the bins (with labels consistent with the physical data)
appsT$dateTimeB <- cut(appsT$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# compute total abundance per taxon per bin
appsB <- ddply(appsT, ~ transect + dateTimeB + group + taxon, function(x) {
  # compute abundance in the bin (actual count * subsampling rate)
  # i.e. as if all organisms were counted on every frame
  abund <- x$count * x$sub
  # compute total abundances
  tot <- sum(abund, na.rm=T)
  # TODO for organisms that were subsampled, this may result in errors because subsampling intervals would overlap with the previous and next frames. This is especially critical for appendicularians which were subsampled every 20 frames while 1m bins correspond to approximately 39 frames when ISIIS is going down fast
  return(c(abund=tot))
}, .progress="text", .parallel=parallel)

# choose only some physical data
phytemp <- phyB[phyB$transect==1 & phyB$cast==2,]
phytemp <- rbind(phytemp, phyB[phyB$transect==1 & phyB$cast==15,])
phytemp <- rbind(phytemp, phyB[phyB$transect==1 & phyB$cast==18,])
phytemp <- rbind(phytemp, phyB[phyB$transect==2 & phyB$cast==18,])
phytemp <- rbind(phytemp, phyB[phyB$transect==2 & phyB$cast==20,])
phytemp <- rbind(phytemp, phyB[phyB$transect==2 & phyB$cast==23,])
phytemp <- rbind(phytemp, phyB[phyB$transect==3 & phyB$cast==3,])
phytemp <- rbind(phytemp, phyB[phyB$transect==3 & phyB$cast==16,])
phytemp <- rbind(phytemp, phyB[phyB$transect==3 & phyB$cast==18,])

phyB <- phytemp

da <- ddply(appsB[,c("dateTimeB", "group", "taxon", "abund")], ~ taxon, function(b, p=phyB){
  # join physical data and counts for this taxon
  da <- join(p, b, type="left", by="dateTimeB")
  # when a taxon is not present in a bin, all columns corresponding to the biological data will contain NA
  # make sure they contain a zero abundance for this taxon instead
  d$group <- b$group[1]
  d$taxon <- b$taxon[1]
  d$abund[is.na(d$abund)] <- 0
  return(d)
}, .progress="text", .parallel=parallel)

# compute concentrations
da$concentration <- da$abund / da$volume

da <- da[da$down.up == "down",]

ggplot(da) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70","black")) + scale_area(range=c(1,10))

# save as csv
write.csv(da, "data/apps_binned_by_depth.csv", row.names=FALSE)

# }

