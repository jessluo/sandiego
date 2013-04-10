#
#      Read and preprocess ISIIS particle count data
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())
parallel <- TRUE

`%ni%` <- Negate(`%in%`) 

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read in file list and cull list of empty files ------------------------

# list files
pcFiles1 <- list.files("raw_particlecount_data/101510", full=TRUE)
pcFiles2 <- list.files("raw_particlecount_data/101610", full=TRUE)

# filter out the files that are smaller or equal to 3 bytes in size (empty files)
pcFiles1 <- pcFiles1[file.info(pcFiles1)$size>3]
pcFiles2 <- pcFiles2[file.info(pcFiles2)$size>3]

# }

##{ Import the data -------------------------------------------------------

# read in files
pc1 <- adply(pcFiles1, 1, function(file){
  
  d <- read.delim(file, header=TRUE, stringsAsFactors=FALSE)
  return(d)

}, .progress="text")

pc2 <- adply(pcFiles2, 1, function(file){
  
  d <- read.delim(file, header=TRUE, stringsAsFactors=FALSE)
  return(d)
  
}, .progress="text")

# remove adply column
pc1 <- pc1[,-1]
pc2 <- pc2[,-1]

pc <- rbind(pc1, pc2)

# save it as a csv file
write.csv(pc, file="raw_particlecount_data/raw_import.csv", row.names=FALSE)

# }

##{ Calculate and create a proper time format ----------------------------

# import particle count data, if not already done so
pc <- read.csv("raw_particlecount_data/raw_import.csv", header=TRUE, stringsAsFactors=FALSE)

# define a beginning position for subsetting the name
# this cuts out the name of the file folder as well as "Dataset3_"
BEGSUB <- 10

# original name from the file
pc$dateName <- str_sub(pc$Label, BEGSUB, str_length(pc$Label)-4)

# calculating the time stamp of the stack
# subsetting just the timestamp
pc$origStack <- str_sub(pc$dateName, 1, 17)
# adding in a period to help R recognize that the seconds are in decimals
pc$origStack <- str_c(str_sub(pc$origStack, 1,14), ".", str_sub(pc$origStack, 15,17))
# striptime and convert to POSIXct format
pc$origStack <- as.POSIXct(strptime(pc$origStack, format="%Y%m%d%H%M%OS", tz="America/Los_Angeles"))

# add the fraction of the stack
pc$frameFr <- str_sub(pc$dateName, 23, str_length(pc$dateName))
pc$frameFr <- as.numeric(pc$frameFr)
pc$dateTime <- pc$origStack + 17.55 * pc$frameFr / 300

# shift time back 3 hours because we recorded it in eastern time
pc$dateTime <- pc$dateTime - 3 * 3600

pc <- pc[, c("X", "dateTime", "Area", "Perim.", "BX", "BY", "Width", "Height", "Major", "Minor", "Angle")]

write.csv(pc, "raw_particlecount_data/particle_count_corr_timestamp.csv", row.names=FALSE)

# }
