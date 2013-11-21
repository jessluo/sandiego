#
#      Bin data (both physical and biological) over depth intervals
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("ggplot2")
# library("oce")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- TRUE

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards
bio <- read.csv("data/bio.csv", stringsAsFactors=FALSE)
bio$dateTime <- as.POSIXct(bio$dateTime, tz="GMT")


##{ Compute time bins on physical data ------------------------------------

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
timeBins <- as.POSIXct(timeBins, origin="1970-01-01 00:00:00.00 GMT", tz="GMT")

# verify where the cuts are
dd <- data.frame(dateTime=timeBins)
# find the depth/long/transect number where the cuts will happen
dd$depth <- approx(phy$dateTime, phy$depth, timeBins)$y
dd$long <- approx(phy$dateTime, phy$long, timeBins)$y
dd$transect <- approx(phy$dateTime, phy$transect, timeBins)$y
# plot the top and bottom of the yos
p <- ggplot(mapping=aes(x=long, y=-depth)) + geom_point(data=phy[phy$transect==1,], size=0.5) + facet_grid(transect~.) + geom_point(data=dd[dd$transect==1,], size=2)
p + coord_cartesian(ylim=c(-4, 0))
p + coord_cartesian(ylim=c(-138, -120))

# }


##{ Bin physical data -----------------------------------------------------

# create the bins (and ensure labels are simple, consistent and unique)
phy$dateTimeB <- cut(phy$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# convert from factor to numeric
phy$dateTimeB <- as.numeric(phy$dateTimeB)

# calculate the time duration of each bin
# assign the timeBin values to the phy dataframe, convert to POSIX
phy$timeBins <- cut(phy$dateTime, breaks=timeBins, labels=timeBins[1:length(timeBins)-1])
phy$timeBins <- as.POSIXct(phy$timeBins, origin="1970-01-01 00:00:00.00 GMT", tz="GMT")

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


# change the NAs - the first set of NAs is from the beginning of the transect to the first bin, and the second set of NAs are from the last bin to the end
# the first set of NAs will be changed to be bin 0 and the last set of NAs will be grouped with the previous bin

# call the first datetimeB "0" instead of NA
# NB this is technically not how you're supposed to do this subsetting/resetting of values but since these values are at the beginning of the dataframe it should be ok
phy[which(is.na(phy[phy$transect==1 & phy$cast==1,"dateTimeB"])),"dateTimeB"] <- 0

phy[which(is.na(phy$dateTimeB)),"dateTimeB"] <- as.numeric(max(phy$dateTimeB, na.rm=T))

# add another dateTimeB column
phy$dateTimeB2 <- phy$dateTimeB

# compute average properties per bin
# TODO because we compute things by date bins but *also* cast and down.up, we are actually separating the descending and ascending part of the turn at depth, for example. There is no reason to think they would be different. We need the separation here because we only have downcasts in the biological data anyway but with the complete data, the binning should be by date bin only.
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
  # try to compute it from distance travelled
  # dist <- max(geodDist(x$lat, x$long, alongPath=TRUE)) * 1000 * 1000
  # (h.vel <- dist / out$duration)
  # # the problem is that there is not enough precision in the lat/long in order to truly calculate a distance
  # assign the mean/median velocity instead
  # choose the mean/median velocity of downcasts
  # median = 2511 mm/s (2.551 m/s)
  # mean = 2511.884 mm/s (2.511884 m/s)
  if (is.na(out$velocity)) {
    out$velocity <- 2.5119
  }
 
  # compute volume sampled in that bin (in m^3)
  #             ISIIS field in m^2  *  speed in m/s           *  duration in s  
  out$volume <- 0.13 * 0.50         *  (out$velocity)         *  out$duration
  # To consider doing: this theoretical field of view is for large organisms like jellies. For smaller appendicularians, your field of view would be smaller (possibly 13 cm x 40 or 45 cm) because of the difficulty in resolving appendicularians when they are too close to the glass.
  
  # TODO check compared to computation using volume per frame and frames per second
  # --> checked with calculating of frames per second, some of the data bins have a duration per frame value greater than 0.75 seconds. using 1 m depth bins, some of the bins have values greater than 1.5 seconds, which seem high. using 2 m depth bin, this smooths out the high values a bit (down to 1.3 seconds per frame) However, this is only present in 1-1.3% of the bins and this percentage seems to be consistent with depth bin sizes up to 5 m
  
  return(out)
}, .progress="text", .parallel=parallel)

# compute total sampling volume
sum(phyB[phyB$transect==1,]$volume)
# [1] 4001.259

sum(phyB[phyB$transect==2 & phyB$cast >=5,]$volume)
# [1] 4427.926

sum(phyB[phyB$transect==3,]$volume)
# [1] 2696.185


# }


##{ Bin biological data ---------------------------------------------------

# create the bins (with labels consistent with the physical data)
bio$dateTimeB <- cut(bio$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# TODO there is an issue with the timeBins in the biological data. NA values in downcast 1, transect 1
bio[which(is.na(bio$dateTimeB)),]
# --> these NAs are data points near the surface where phy data have NA values in dateTimeB (first depth bin is not marked)

# until can fix the binning issue, remove these data points
bio <- bio[complete.cases(bio),]

# compute total abundance per taxon per bin
bioB <- ddply(bio, ~ transect + cast.bio + dateTimeB + group + taxon, function(x) {
  # compute abundance in the bin (actual count * subsampling rate)
  # i.e. as if all organisms were counted on every frame
  abund <- x$count * x$sub
  # compute total abundances
  tot <- sum(abund, na.rm=T)
  # TODO for organisms that were subsampled, this may result in errors because subsampling intervals would overlap with the previous and next frames. This is especially critical for appendicularians which were subsampled every 20 frames while 1m bins correspond to approximately 39 frames when ISIIS is going down fast
  return(c(abund=tot))
}, .progress="text", .parallel=parallel)

# }


##{ Join physical and biological data -------------------------------------
# take physical data as reference to create zeros when nothing was caught in a bin

# ---this section is commented out because the selection shouldn't be done manually before the binning---
# select only downcasts since biological data was aquired only on down casts
# phyB <- phyB[phyB$down.up == "down",]

# one problem with selecting only downcasts is that if one depth bin spans both downcasts and upcasts, then both the duration of the depth bin as well as the volume sampled is wrong, which changes the calculated concentration of the organisms

# join bio and phy data by time bin
# make sure we create zeros for species absent in each bin
d <- ddply(bioB[,c("dateTimeB", "group", "taxon", "abund")], ~ taxon, function(b, p=phyB){
  # join physical data and counts for this taxon
  d <- join(p, b, type="left", by="dateTimeB")
  # when a taxon is not present in a bin, all columns corresponding to the biological data will contain NA
  # make sure they contain a zero abundance for this taxon instead
  d$group <- b$group[1]
  d$taxon <- b$taxon[1]
  d$abund[is.na(d$abund)] <- 0
  return(d)
}, .progress="text", .parallel=parallel)

# compute concentrations
d$concentration <- d$abund / d$volume

# TODO: remove this later when you have solved the NA problem in computing timeBins
d <- d[complete.cases(d$dateTimeB),]


# remove the zeros for the transects / casts in which data was not recorded (these are not true zeros)
# steps:
# 1. Use ddply to mark and delete transects in which no data was recorded for each group. 
d <- ddply(d, ~ transect + group, function(x){
  x$tfzero <- sum(x$abund) != 0
  return(x)
})

d <- d[which(d$tfzero),]
# problem here because some hydromedusae and appendicularians somehow got merged with physical data from the third transect, which should not happen (hydromedusae and appendicularians were not counted yet in the 3rd transect) 
# this has to do with the NAs in the timeBins ... see above comment in the timeBins section. Maybe caused by the binning process?

# 2. Mark the depth bins in which both down casts and upcasts are captured
# this can't be done here --- must be done above before binning

# 3. Remove all upcasts 
d <- d[which(d$down.up == "down"),]

# 4. Use ddply to identify and mark all down casts in which there is no data for each group
d <- ddply(d, ~transect + cast + group, function(x){
  x$tfzero <- sum(x$abund) != 0
  return(x)
})

d <- d[which(d$tfzero),]

# inspect
d1 <- d[which(d$tfzero==F),]

d_ply(d1, ~transect + cast + group, function(x){
  s <- c(unique(x$transect), unique(x$cast), unique(x$group))
  print(s)
})
 
d <- d[,-which(names(d)=="tfzero")]# 
# }


##{ Inspect joined data ---------------------------------------------------

summary(d)

ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=abund)) + facet_wrap(~group, scales="free")
ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group, scales="free")
ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=abund)) + facet_wrap(~taxon, scales="free")
ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=concentration)) + facet_wrap(~taxon, scales="free")

ggplot(d) + geom_histogram(aes(x=volume), binwidth=0.1)

# inspect concentrations and print them out into a PDF
# pdf("binned_concentrations.pdf", width=14, height=10)
# 
alply(unique(d$group), 1, function(group) {
	g <- ggplot(d[d$group==group,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + 
  facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_size_area(max_size=12)
  print(g)
})
# dev.off()

d$front <- factor(d$front)
d$front2 <- d$front
d$front2[d$front2 !="front"] <- "not_front"
d$front2 <- factor(d$front2)

for(i in 1:length(unique(d$taxon))) {
   taxon <- unique(d$taxon)[i]
    print(taxon)
    print(kruskal.test(concentration~front, data=d[d$taxon==taxon & d$front !="west",]))
}

# print a the range of concentrations into a PDF
library(gridExtra)
pdf("orig_concentration_range.pdf", width=7, height=14)
print_conc <- adply(unique(d$taxon), 1, function(taxon)
  {
      df <- data.frame(taxon, min=min(d[d$taxon==taxon,"concentration"]), max=max(d[d$taxon==taxon,"concentration"]))
      return(df)
  })
print_conc <- print_conc[,-1]
grid.table(print_conc)
dev.off()

# inspect abundances
alply(unique(d$group), 1, function(group) {
  ggplot(d[d$group==group,]) + geom_point(aes(x=long, y=-depth, size=abund, colour=abund>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))
})

pdf("group_abund.pdf", width=8.5, height=11)
print_conc <- adply(unique(d$group), 1, function(group)
{
  df <- data.frame(group, conc=mean(d[d$group==group,"concentration"]), abund=sum(d[d$group==group,"abund"]))
  return(df)
})
print_conc <- print_conc[,-1]
grid.table(print_conc)
dev.off()

# subset plots
ggplot(d[d$taxon %in% c("h5_Liriope", "h6_Solmundella", "vsh", "h7_Pegantha", "h7_Rhopalonema", "h9_Aglaura", "h2_Haliscera"),]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

ggplot(d[d$taxon=="lemu",]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

ggplot(d[d$taxon %in% c("Beroida", "Hormiphora californiensis", "Larval Lobata", "Ocyropsis maculata", "Velamen"),]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=concentration>0), alpha=0.7) + facet_grid(transect~taxon) + scale_colour_manual(values=c("grey70", "black")) + scale_area(range=c(1,10))

ggplot(d[d$taxon %in% c("Larval Lobata", "Ocyropsis maculata") & d$concentration > 0,]) + geom_point(aes(x=long, y=-depth, size=concentration, colour=taxon), alpha=0.7) + facet_grid(transect~.) + scale_area(range=c(1,10))
# }

write.csv(d, "data/all_binned_by_depth.csv", row.names=FALSE)
