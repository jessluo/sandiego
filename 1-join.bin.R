#
#   Joining and Binning ISIIS San Diego data
#
#   Jean-Olivier Irisson and Jessica Luo
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Joining physical and biological data ----------------------------------

# start from scratch just to be sure
bio <- read.csv("data/bio.csv", colClasses=c(dateTime="POSIXct"), stringsAsFactors=FALSE)
phy <- read.csv("data/phy.csv", colClasses=c(dateTime="POSIXct"), stringsAsFactors=FALSE)

# subset the physical data to match where the biological data was recorded
# TODO generalize the subseting of data (actually there should not be any subseting at all in the end ;) )
# which biological downcasts have all the data?
dc <- unique(bio$cast.bio[which(bio$group=="Hydromedusae" & bio$transect==2)])
# select those downcasts
phy <- phy[which(phy$down.up=="down" & phy$transect==2 & phy$cast %in% dc),]
bio <- bio[which(bio$transect==2 & bio$cast.bio %in% dc),]

# Bin bio and physical data by time, on ~ 20 frames
# appendicularians are looked for every 20 frames and we want all bins to (potentially) have appendicularian counts
# once this binning is done we can compute true abundances (count * subsampling rate) and finally bin over larger scale bins

bio$dateTimeB1 <- round_any(bio$dateTime, 1.2)
phy$dateTimeB1 <- round_any(phy$dateTime, 1.2)

# compute total abundance in each time bin
bioB1 <- ddply(bio, ~ transect + cast.bio + dateTimeB1 + group + taxon + sub, function(x){
  return(c(count=sum(x$count)))
}, .progress="text")
# scale that count by the subsetting interval
bioB1$abund <- bioB1$count * bioB1$sub

# compute average/total of physical properties in each time bin
phyB1 <- ddply(phy, ~ transect + cast + down.up + dateTimeB1, function(x){
  # means of columns that are numbers
  means <- colMeans(x[,llply(phy, class) == "numeric"])
  # TODO should add sum of vol.imaged when that works
  # compute time elapsed 
  return(means)
}, .progress="text")

# join bio and phy data by time, take the physical data as a reference
# make sure we create zeros for species absent in each bin
d <- ddply(bioB1[,c("dateTimeB1", "group", "taxon", "abund")], ~ taxon, function(b, p=phyB1){
  # join physical data and counts for this taxon
  d <- join(p, b, type="left", by="dateTimeB1")
  # when a taxon is not present in a bin, all columns corresponding to the biological data will contain NA
  # make sure they contain a zero abundance for this taxon instead
  d$group <- b$group[1]
  d$taxon <- b$taxon[1]
  d$abund[is.na(d$abund)] <- 0
  return(d)
}, .progress="text")

# reorder the data by time, just to be cleaner
d <- arrange(d, dateTimeB1, group, taxon)

# write that data set
write.csv(d, "data/full_bin1.csv", row.names=FALSE)

# }

##{ Binning the data by depth ---------------------------------------------

d <- read.csv("data/full_bin1", colClasses=c(dateTimeB1="POSIXct"), stringsAsFactors=FALSE)

# what is a reasonable depth to bin by? 1 meter?
# 1 meter depth change traveling down at 10 degrees ~ 5.75 of distance covered, and traveling at 2.5 m/s that is 2.304 sec. which is approx 39 images
# traveling at 5 degrees, distances are doubled, so integrating over 78 images
# 1 meter depth bins are probably the minimum

d$depthB2 <- round_any(d$depth, 1)

# TODO by ddplying on cast + depthB2 we are actually separating the descending and ascending part of the turn at depth, for example, while there is no reason to think they would be different. This is not a problem here since we only have downcasts but in the complete data, the binning should be by depth only (and therefore the computation cannot be done with just ddply)
dB2 <- ddply(d, ~ cast + depthB2 + taxon, function(x){
  # total abundance
  abund <- sum(x$abund)
  # location
  loc <- colMeans(x[,c("lat", "long")])
  # time
  dateTime <- mean(x$dateTimeB1)
  # physical data
  phy <- colMeans(x[,c("temp", "salinity", "fluoro", "oxygen", "irrandiance", "heading")])
  
  # count number of "seconds" (that's actually the nb of 1 sec bins)
  # TODO we should be more precise here...
  duration <- nrow(x)
  
  return(data.frame(dateTime, t(loc), t(phy), abund, duration))
}, .progress="text")

# so assuming 17 frames per second and field of view of 13 cm x 45 cm, ISIIS images 0.1293 m^3/s. 
# so equivalent time to image 1 m^3 is 7.735 sec. equivalent to 131.5 frames.
# physical and biological data is already binned at 1 sec intervals. 
# so if you divide by binBy, then you will get counts per 0.1293 m^3. So just multiply by the inverse and get counts per m^3.
# TODO not sure I understand the reasonning here. It must depend on the speed of trawling...
dB2$density <- dB2$count / dB2$duration * 7.735

write.csv(dB2, "data/full_bin2.csv", row.names=FALSE)

# }
