#
#      Read and preprocess ISIIS particle count data
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("ggplot2")
# library("foreach")
# library("doParallel")
# registerDoParallel(cores=detectCores())
# parallel <- TRUE

`%ni%` <- Negate(`%in%`) 

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read in folder list ---------------------------------------------------

# list folders
pcFolders1 <- list.files("../SD_ImageJ/particleCount_2013-03/101510", full=TRUE)
pcFolders2 <- list.files("../SD_ImageJ/particleCount_2013-03/101610", full=TRUE)


# }

##{ Import the data -------------------------------------------------------

d1 <- data.frame()
d2 <- data.frame()

for (i in 1:length(pcFolders1)) {
  
  # initialize
  
  # get a list of files within each folder
  pcFiles <- list.files(pcFolders1[i], full=TRUE)
  
  # filter out the files that are smaller than 3 bytes in size (empty files)
  pcFiles <- pcFiles[file.info(pcFiles)$size >3]
  
  # read in the files
  for (j in 1:length(pcFiles)){
    di <- read.delim(pcFiles[j], header=TRUE, stringsAsFactors=FALSE)
    
    # combine the files
    d1 <- rbind(d1, di)
  }
  print(i)
}

for (i in 1:length(pcFolders2)) {
  
  # get a list of files within each folder
  pcFiles <- list.files(pcFolders2[i], full=TRUE)
  
  # filter out the files that are smaller than 3 bytes in size (empty files)
  pcFiles <- pcFiles[file.info(pcFiles)$size >3]
  
  # read in the files
  for (j in 1:length(pcFiles)){
    di <- read.delim(pcFiles[j], header=TRUE, stringsAsFactors=FALSE)
    
    # combine the files
    d2 <- rbind(d2, di)
  }
  print(i)
}

pc <- rbind(d1, d2)

# save it as a csv file
write.csv(pc, file="data/raw_particle_count_import.csv", row.names=FALSE)

# }

##{ Calculate and create a proper time format ----------------------------

# import particle count data, if not already done so
pc <- read.csv("data/raw_particle_count_import.csv", header=TRUE, stringsAsFactors=FALSE)

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

write.csv(pc, "data/particle_count_corr_timestamp.csv", row.names=FALSE)

# }
