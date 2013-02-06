#
#   Read and preprocess ISIIS data
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

bio <- read.csv("data/bio.csv", colClasses=c(dateTime="POSIXct"), stringsAsFactors=FALSE)
phy <- read.csv("data/phy.csv", colClasses=c(dateTime="POSIXct"), stringsAsFactors=FALSE)

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
