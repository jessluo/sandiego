#
#      Read, process and add physical data for appendicularian data
#      Collected by Jenna Binstein
#         Data collected on transects 1-2 of San Diego Data
#         Downcasts: 5, 6, 12, 13, 16, 21 of transect 1
#                    5, 11, 25, 27, 29, 30 of transect 2
#         Apps counted in 1/5 of every 20th frame
#
#      Code by Jessica Luo, 2013
#
#--------------------------------------------------------------------------

library(plyr)
library(reshape2)
library(ggplot2)
library(stringr)

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read in physical and biological data ----------------------------------
# read physical data
phy <- read.csv("data/phy.csv", stringsAsFactors=FALSE)
phy$dateTime <- as.POSIXct(phy$dateTime, tz="America/Los_Angeles")

# read in appendicularian data
app <- read.csv("raw_biological_data/Appendicularians_JennaBinstein.csv", na.strings="-", stringsAsFactors=FALSE)

# clean up the names
names(app) <- tolower(names(app))
app <- rename(app, c("lowest.classification" = "taxon", "inshore.offshore.front" = "location"))

# TODO: enumerate the transects instead of having it be specified by date

# }

##{ Calculate and create a proper time format ----------------------------

app$year <- str_c("20", str_sub(app$transect, 7, 8))
app$month <- str_sub(app$transect, 1, 2)

# date moves ahead one unit when the time crosses zero
app$date <- ifelse(app$hour >= 18 & app$hour <= 23, str_sub(app$transect, 4, 5), as.numeric(str_sub(app$transect, 4, 5)) + 1)

app$sec <- str_c(app$sec, ".", app$s.1000)
# recompute the time
# the record contains the time of the first frame in the stack and the frame number within the stack
# there are 300 frames in a stack and they represent a time of 17.55 seconds
app$dateTime <- as.POSIXct(str_c(app$year, "-", app$month, "-", app$date, " ", sprintf("%02i",app$hour), ":", sprintf("%02i",app$min), ":", app$sec), tz="America/Los_Angeles")
app$dateTime <- app$dateTime + 17.55 * app$frame / 300

# shift time back 3 hours because we recorded it in eastern time
app$dateTime <- app$dateTime - 3 * 3600

# keep relevant columns
app <- app[, c("transect", "downcast", "location", "file", "dateTime", "taxon", "in.out", "size", "degradation")]

# rename transect column
app <- rename(app, c("transect"="transectDate"))
# }


##{ Interpolate physical data ---------------------------------------------
# interpolate physical parameters for each biological sample

# unique biological sample times
times <- sort(unique(app$dateTime))

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

# put all that back in img
app <- join(app, interp)

# }

##{ Keep columns and save as csv ------------------------------------------

# columns to keep
app <- app[,c("transectDate", "downcast", "location", "file", "dateTime", "taxon", "in.out", "size", "degradation", "depth", "lat", "long", "temp", "salinity", "fluoro", "oxygen")]

app <- rename(app, c("transectDate" = "transect"))

# write it as text file
write.csv(app, file="data/jbinstein_apps_with_phys_vars.csv", row.names=FALSE)

# }
