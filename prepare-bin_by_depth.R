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
library("PBSmapping")
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

# per cast, interpolate the times corresponding to the crossing of the depth bins
timeBins <- dlply(phy, ~transect+cast+down.up, function(x, bin) {
  # maximum depth reached in this cast
  maxD <- max(x$depth, na.rm=T)
  
  # compute depth bins crosses
  depthBins <- seq(0, maxD, by=bin)
  
  # compute corresponding crossing times
  timeBins <- approx(x$depth, x$dateTime, depthBins)$y
  # TODO try to extrapolate to get the first points and not loose data
  
  return(timeBins)
}, bin=binSize, .progress="text")
# combine the elements of the lists into just one vector
timeBins <- sort(do.call(c, timeBins))

# make it a proper time
timeBins <- as.POSIXct(timeBins, origin="1970-01-01 00:00:00.00 GMT", tz="GMT")

# }


##{ Bin physical data -----------------------------------------------------

# create the bins (and ensure labels are simple, consistent and unique)
phy$dateTimeB <- cut(phy$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# compute average properties per bin
# TODO because we compute things by date bins but *also* cast and down.up, we are actually separating the descending and ascending part of the turn at depth, for example. There is no reason to think they would be different. We need the separation here because we only have downcasts in the biological data anyway but with the complete data, the binning should be by date bin only.
phyB <- ddply(phy, ~transect+cast+down.up+dateTimeB, function(x) {
  # means of variables inside the bin
  out <- as.data.frame(as.list(colMeans(x[,c("depth", "lat", "long", "temp", "salinity", "fluoro", "oxygen", "irradiance", "horizontal.vel", "vertical.vel", "pitch", "velocity")], na.rm=T)))
  
  # compute duration (in seconds)
  out$duration <- as.numeric(diff(range(x$dateTime, na.rm=T)))
  
  # recompute vertical velocity in this larger bin (more precise)
  # remove records in which either depth or time is missing
  xS <- x[,c("depth", "dateTime")]
  xS <- na.omit(xS)
  # compute the vertical speed
  out$vertical.vel <- diff(range(xS$depth)) / as.numeric(diff(range(xS$dateTime))) * 1000

  # TODO recompute horizontal velocity from the lat/long difference when it is NA/NaN
  #      check out geodDist in oce to compute distances from lat and long
  dist <- geodDist(phy$lat, phy$long, alongPath=TRUE)
  # computes the distance differences along the path
  dist <- diff(dist)
  # the problem is that there is not enough precision in the lat/long in order to truly calculate a distance
  
  # you could also calculate a euclidean distance between UTM coordinates + depth (and then obtain the total distance between points)
  xy<-cbind(X=phy$long, Y=phy$lat) # longitude is in the x direction, latitude is in the y direction
  attr(xy,"zone")=11 # zone where this study site is located
  attr(xy, "projection")="LL"
  # convert to UTM
  XY<-convUL(xy)
  UTM.X<-XY[,1]
  UTM.Y<-XY[,2]
  # does plain "dist" work? apparently not, there is an error about negative length vectors ...
  # dist2 <- dist(UTM.X, UTM.Y, -phy$depth, method="euclidean") 
  
  
  # recompute velocity since we recomputed vertical velocity
  out$velocity <- sqrt(out$horizontal.vel^2 + out$vertical.vel^2)

  # compute volume sampled in that bin (in m^3)
  #             ISIIS field in m^2  *  speed in m/s           *  duration in s  
  out$volume <- 0.13 * 0.45         *  (out$velocity / 1000)  *  out$duration
  # TODO check compared to computation using volume per frame and frames per second

  return(out)
}, .progress="text", .parallel=parallel)

# }


##{ Bin biological data ---------------------------------------------------

# create the bins (with labels consistent with the physical data)
bio$dateTimeB <- cut(bio$dateTime, breaks=timeBins, labels=1:(length(timeBins)-1))

# compute total abundance per taxon per bin
bioB <- ddply(bio, ~ transect + cast.bio + dateTimeB + group + taxon, function(x) {
  # compute abundance in the bin (actual count * subsampling rate)
  # i.e. as if all organisms were counted on every frame
  abund <- x$count * x$sub
  # compute total abundances
  tot <- sum(abund, na.rm=T)
  # TODO for organisms that were subsampled, this may result in errors because subsampling intervals would overlap with the previous and next frames. This is especially critical for appendicularians which were subsampled every 20 frames while 1m bins correspond to approximately 39 frames when ISIIS is goind down fast
  return(c(abund=tot))
}, .progress="text", .parallel=parallel)

# }


##{ Join physical and biological data -------------------------------------
# take physical data as reference to create zeros when nothing was caught in a bin

# select only downcasts since biological data was aquired only on down casts
# phyB <- phyB[phyB$down.up == "down",]

# select only downcasts in which all the biological data is present (transect 2 only)
dc <- unique(bio$cast.bio[which(bio$group=="Hydromedusae" & bio$transect==2)])
phyB <- phyB[which(phyB$down.up=="down" & phyB$transect==2 & phyB$cast %in% dc),]
bioB <- bioB[which(bioB$transect==2 & bioB$cast.bio %in% dc),]

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

# }


##{ Inspect joined data ---------------------------------------------------

summary(d)

ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=abund)) + facet_wrap(~group, scales="free")
ggplot(d[d$abund != 0,]) + geom_histogram(aes(x=abund)) + facet_wrap(~taxon, scales="free")

for (group in unique(d$group)) {
  print(ggplot(d[d$group==group,]) + geom_point(aes(x=long, y=-depth, size=abund, colour=abund>0), alpha=0.7) + facet_grid(taxon~transect) + scale_colour_manual(values=c("grey70", "black")))
}
# -> there seem to be a few very high values of abundance, often near the surface
# TODO find out wether this is present in the original data or is caused by the binning

# }

write.csv(d, "data/all_binned_by_depth.csv", row.names=FALSE)
