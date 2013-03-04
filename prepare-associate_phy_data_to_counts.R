#
#      Associate physical variables to each count through interpolation
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

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

# read data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards
bio <- read.csv("data/bio.csv", stringsAsFactors=FALSE)
bio$dateTime <- as.POSIXct(bio$dateTime, tz="GMT")


##{ Interpolate physical data ---------------------------------------------
# interpolate physical parameters for each biological sample

# unique biological sample times
times <- sort(unique(bio$dateTime))

# prepare storage
interp <- data.frame(dateTime=times)

# interpolate all numeric columns linearly
cols <- names(phy)[which(llply(phy, class)=="numeric")]
for (col in cols) {
  interp[,col] <- approx(phy$dateTime, phy[,col], times)$y
}

# get transect, cast etc. using nearest neighbour interpolation
closestIndex <- round(approx(phy$dateTime, 1:nrow(phy), times)$y)
interp <- cbind(interp, phy[closestIndex, c("transect", "cast", "down.up")])

# put all that back in bio
bio <- join(bio, interp)

# }

# convert counts into presence
bio$presence <- bio$count > 0

# write it as text
write.csv(bio, file="data/all_presence_only.csv", row.names=FALSE)
