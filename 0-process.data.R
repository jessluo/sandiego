#
#   Read and preprocess ISIIS data
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)


##{ Physical data ---------------------------------------------------------

# list all the physical data files in a given directory
phyFiles <- list.files("raw_physical_data", full=TRUE)

# read them all
phy <- adply(phyFiles, 1, function(file) {
  
  # read the data
  d <- read.table(file, sep="\t", skip=10, header=TRUE, fileEncoding="ISO-8859-1", stringsAsFactors=FALSE, quote="\"", check.names=FALSE, encoding="UTF-8", na.strings="9999.99")
  
  
  # clean names
  head <- names(d)
  head <- str_replace(head, "\\(.*\\)", "")
  head <- str_trim(head)
  head <- make.names(head)
  head <- tolower(head)
  head <- str_replace(head, fixed(".."), ".")
  # assign names
  names(d) <- head
  
  # create a proper date + time format
  date <- scan(file, what="character", skip=1, nlines=1, quiet=TRUE)
  date <- date[2]
  mm <- str_sub(date,1,2)
  dd <- str_sub(date,4,5)
  dd <- as.numeric(dd)
  yy <- str_sub(date,7,8)
  dateNextDay <- str_c(mm,as.character(dd+1),yy, sep="/")
  
  # shift by one day when we cross midnight
  d$hour <- as.numeric(str_sub(d$time,1,2))
  d$date <- ifelse(d$hour >= 18 & d$hour <= 23, date, dateNextDay)
  d$dateTime <- str_c(d$date, d$time, sep=" ")
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="UTC"))
  
  # shift all physical data back 3 hours
  d$dateTime <- d$dateTime - 3 * 3600
  
  # code in a transect number
  # this is not robust for all physical data but is necessary here
  d$transect <- dd-14
  
  # reformat the lat and long in decimal degrees
  to.dec <- function(x) {
    # split in degree, minute, second
    pieces <- str_split_fixed(x, "° |'", 3)
    # extract orientation (S/N and E/W)
    orientation <- str_sub(pieces[,3], -1)
    # remove orientation to only keep numbers
    pieces[,3] <- str_replace(pieces[,3], "[SNWE]", "")
    # convert to decimal degrees
    dec <- as.numeric(pieces[,1]) + as.numeric(pieces[,2]) / 60 + as.numeric(pieces[,3]) / 3600
    # orient the coordinate
    ifelse(orientation %in% c("S", "W"), -dec, dec)
    
    return(dec)
  }
  d$lat <- to.dec(d$lat)
  d$long <- to.dec(d$long)
  
  # keep only interesting data
  d <- d[,c("dateTime", "depth", "lat", "long", "temp", "salinity", "fluoro", "oxygen", "irrandiance", "vol.imaged", "heading", "transect")]
  
  # write it as a CSV file
  # outFile <- basename(str_replace(file, "txt", "csv"))
  # outFile <- str_c("data/", outFile)
  # write.csv(d, file=outFile, row.names=FALSE) 
  
  return(d)
  
}, .progress="text")

# remove adply crap
phy <- phy[,-1]

# Detect up and down casts and number them
phy <- ddply(phy, ~transect, function(d) {
  # detect up and down casts by smoothing the depth profile and finding the turning points
  # smooth depths
  order <- 100
  # tested options between 10-20, correct number after order=14. use order between 16-20 for robustness
  depth_avg <- decaverage(-d$depth, times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
  # plot(depth_avg)
  depth_avg <- as.numeric(extract(depth_avg, component="filtered"))
  # detect turning points
  TP <- suppressWarnings(turnpoints(depth_avg))
  
  # NB: if we were assigning a different number of every up and down cast we would do it like that
  # castNb <- cumsum(TP$peaks | TP$pits)
  # castType <- castNb %% 2
  # but we don't do this here
  
  # number casts; a cast is everything between two turns (either at the surface or at depth)
  # detect wether the first turn is at the surface or at depth
  turnAtSurface <- which(TP$peaks)
  turnAtDepth <- which(TP$pits)
  firstTurn <- c("surface", "depth")[which.min(c(turnAtSurface[1], turnAtDepth[1]))]
  
  # separate up and down casts
  castAbsNb <- cumsum(TP$peaks | TP$pits)
  castType <- castAbsNb %% 2 + 1
  
  # assign cast numbers and types based on the first turn
  if (firstTurn == "surface") {
    # when the first turn is at the surface, we started at depths and we warn to number our casts between turns at depth (i.e. pits)
    d$cast <- cumsum(TP$pits) + 1
    # and we start with an upcast
    d$down.up <- c("up", "down")[castType]
  } else {
    # and the contrary when the first turn is at depth
    # NB: this is the most common scenario
    d$cast <- cumsum(TP$peaks) + 1     
    d$down.up <- c("down", "up")[castType]
  }
  
  return(d) 
}, .progress="text")

# visualize
# ggplot(data=phy) + geom_path(aes(x=long, y=-depth, colour=as.factor(cast), linetype=factor(down.up))) + facet_grid(transect ~.) 

# save it as text
write.csv(phy, "data/phy.csv", row.names=FALSE)

# }


##{ Biological data -------------------------------------------------------

bioData <- "raw_biological_data"

# Solmaris
solFiles <- list.files(bioData, pattern="Solmaris", full=TRUE)
sol <- adply(solFiles, 1, read.csv, stringsAsFactors=TRUE)

# identify transects (the index corresponding to the different files)
sol$X1 <- solFiles[sol$X1]
sol$X1 <- str_replace(sol$X1, bioData, "")
sol$transect <- as.numeric(str_sub(sol$X1, 4, 5)) - 14

# create a true date+time column
sol$dateTime <- as.POSIXct(str_c(sol$date, " ", sprintf("%02i",sol$hour), ":", sprintf("%02i",sol$min), ":", sprintf("%02i",sol$sec), ".", sol$s.1000))

# convert to the tall format
sol <- sol[,!names(sol) %in% c("X1", "date", "hour", "min", "sec", "s.1000")]
solT <- melt(sol, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
solT$taxon <- as.character(solT$taxon)
solT$group <- "Solmaris"
# remove the zeros
solT <- solT[solT$count != 0,]

# add subsampling interval (Solmaris were looked for in one every five frame)
solT$sub <- 5


# Siphonophores
siphFiles <- list.files(bioData, pattern="Siphonophores", full=TRUE)
siph <- adply(siphFiles, 1, read.csv)

# identify transects (the index corresponding to the different files)
siph$X1 <- siphFiles[siph$X1]
siph$X1 <- str_replace(siph$X1, bioData, "")
siph$transect <- as.numeric(str_sub(siph$X1, 4, 5)) - 14

# create a true date+time column
siph$dateTime <- as.POSIXct(str_c(siph$date, " ", sprintf("%02i",siph$hour), ":", sprintf("%02i",siph$min), ":", sprintf("%02i",siph$sec), ".", siph$s.1000))

# combine with and without tail
siph$Type1 <- siph$Type1 + siph$Type1_wotail
siph$Type2 <- siph$Type2 + siph$Type2_wotail
siph$Type3 <- siph$Type3 + siph$Type3_wotail
siph <- siph[,-which(str_detect(names(siph), "_wotail"))]

#give siphonophore taxa meaningful names
siph <- rename(siph, c("Type1"="Diphyidae", "Type2"="Sphaeronectes", "Type3"="Physonect", "Type4"="Prayidae", "Type5"="Lilyopsis"))


# convert to the tall format
siph <- siph[,!names(siph) %in% c("X1", "date", "hour", "min", "sec", "s.1000")]
siphT <- melt(siph, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
siphT$taxon <- as.character(siphT$taxon)
siphT$group <- "Siphonophores"
# remove the zeros
siphT <- siphT[siphT$count != 0,]

# add subsampling interval (Siphonophores were looked for in every frame)
siphT$sub <- 1


# Ctenophores
cteFiles <- list.files(bioData, pattern="Ctenophores", full=TRUE)
cteT <- adply(cteFiles, 1, function(file){
  d <- read.csv(file, stringsAsFactors=FALSE)
  
  # clean up taxon names
  d$species <- str_trim(d$species)
  
  # recompute the time
  # the record contains the time of the first frame in the stack and the frame number within the stack
  # there are 300 frames in a stack and they represent a time of 17.375 seconds
  d$dateTime <- as.POSIXct(str_c(d$year, "-", sprintf("%02i",d$month), "-", sprintf("%02i",d$date), " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec)))
  d$dateTime <- d$dateTime + 17.375 * d$x.300 / 300
  
  # extract downcast number
  downcast <- d$folder
  downcast <- str_sub(downcast, start=4, end=5)
  downcast <- str_trim(downcast)
  downcast <- as.numeric(downcast)
  d$downcast <- downcast
  
  # keep only the columns of interest
  d$group <- "Ctenophores"
  d <- rename(d, c("species"="taxon"))
  d <- d[, c("downcast", "dateTime", "group", "taxon")]
  d$count <- 1
  
  return(d)
})

#correct spelling on one taxon
cteT$taxon[which(cteT$taxon=="Ocyropsis maciluta")] <- "Ocyropsis maculata"
# identify transects (the index corresponding to the different files)
cteT$X1 <- cteFiles[cteT$X1]
cteT$X1 <- str_replace(cteT$X1, bioData, "")
cteT$transect <- as.numeric(str_sub(cteT$X1, 4, 5)) - 14
cteT <- cteT[,-1]

# add subsampling interval (Ctenophores were looked for in every frame)
cteT$sub <- 1


# Hydro + Doliolids + Apps
hFiles <- list.files(bioData, pattern=glob2rx("*Hydromedusae*csv"), full=TRUE)
h <- adply(hFiles, 1, function(file) {
  d <- read.csv(file, stringsAsFactors=FALSE)
  names(d) <- tolower(names(d))
  
  # missing values (blanks) are actually 0
  d[is.na(d)] <- 0
  
  # compute the date + time
  # detect starting day from the file name
  file <- str_replace(file, bioData, "")
  day <- as.numeric(str_sub(file, 4, 5))
  # use that to determine the transect number
  d$transect <- day - 14
  # shift by one day when we cross midnight
  d$day <- ifelse(d$hour >= 18 & d$hour <= 23, day, day+1)
  # convert that into POSIXct
  d$dateTime <- as.POSIXct(str_c("2010-10-", d$day, " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec), ".", d$s.1000))
  # remove extraneous date columns
  d <- d[,!names(d) %in% c("day", "hour", "min", "sec", "s.1000")]
  
  return(d)
})

# remove taxa we are not interested in
h <- h[,!names(h) %in% c("muat", "lemu", "spko", "si.2.eu", "nabi", "agel", "coor")]
# and the stuff that's added by adply
h <- h[,!names(h) %in% c("X1")]

#give hydromedusae meaningful names
h <- rename(h, c("h2"="h2_Haliscera", "h3"="h3_Cunina", "h5a"="h5_Liriope", "h6"="h6_Solmundella", "h7"="h7_Rhopalonema", "h7."="h7_Pegantha", "h9"="h9_Aglaura", "h9."="h9_Arctapodema", "h10"="h10_Pegantha", "h11"="h11_Haliscera","h14"="h14_Pegantha", "annatiara"="Annatiara", "r4.aegina"="r4_Aegina", "r5"="r5_Eutonia"))

# convert to the tall format
hT <- melt(h, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
hT$taxon <- as.character(hT$taxon)
# remove the zeros
hT <- hT[hT$count != 0,]

# extract the appendicularians
appT <- hT[hT$taxon=="appendicularians",]
appT$group <- "Tunicates"

hT <- hT[! hT$taxon %in% c("appendicularians", "app.in.hse"),]
hT$group <- ifelse(hT$taxon == "doliolids", "Tunicates", "Hydromedusae")

# add subsampling interval
# (Hydromedusae were looked for in every frame)
hT$sub <- 1
# (Appendicularians were looked for in 1 fifth of a frame, every 20 frames)
appT$sub <- 100
# (except Appendicularians in downcast 1, looked for in one fifth of a frame, every 10 frames)
appT$sub[which(appT$downcast == 1)] <- 50


# put all the data together
bio <- rbind.fill(solT, siphT, cteT, appT, hT)

# shift the time by 3 hours
bio$dateTime <- bio$dateTime - 3 * 3600

bio <- rename(bio, c("downcast"="cast"))

# check
head(bio)
unique(bio$transect)
unique(bio$cast)
unique(bio$group)
unique(bio$taxon)
range(bio$dateTime)
unique(bio$count)
count(bio$count)

# save it, just in case
write.csv(bio, "bio.csv", row.names=FALSE)

# }

# Question: when reading phy data, how to read dateTime as R time format instead of a string?

##{ Joining physical and biological data ----------------------------------

# bin bio and phy data on ~20 frames (because of app)
# round to the nearest second
bio$dateTimer <- round_any(bio$dateTime, 1)
phy$dateTimer <- round_any(phy$dateTime, 1)

bio <- ddply(bio, ~dateTimer+taxon, function(x){
  countsum <- sum(x$count)
  return(data.frame(x[,names(x) != "count"],"count"=countsum))
})

phy <- ddply(phy, ~dateTimer, function(x){
  depth <- mean(x$depth, na.rm=TRUE)
  lat <- mean(x$lat, na.rm=TRUE)
  long <- mean(x$long, na.rm=TRUE)
  temp <- mean(x$temp, na.rm=TRUE)
  salinity <- mean(x$salinity, na.rm=TRUE)
  fluoro <- mean(x$fluoro, na.rm=TRUE)
  oxygen <- mean(x$oxygen, na.rm=TRUE)
  irradiance <- mean(x$irradiance, na.rm=TRUE)
  vol.imaged <- sum(x$vol.imaged, na.rm=TRUE)
  heading <- mean(x$heading, na.rm=TRUE)
  return(data.frame(x[,names(x) %in% c("dateTimer", "transect", "cast", "down.up")], depth, lat, long, temp, salinity, fluoro, oxygen, irradiance, vol.imaged, heading))
}, .progress="text")


# select some downcasts
phy_all <- phy
phy <- phy[which(phy$down.up=="D" & phy$transect==2),]
# which downcasts have all the data?
dc <- unique(bio$cast[which(bio$group=="Hydromedusae" & bio$transect==2)])
dc <- sort(dc, FALSE) #these are the groups we are using
phy <- phy[which(phy$cast %in% dc),]
#yes technically this works BUT the downcasts are ordered wrong because the physical data was broken apart

# join bio and phy data by time

data <- join(phy, bio)
#take the physical data as a reference and fill the missing counts with 0 (see join options)

# }


##{ Binning the data ------------------------------------------------------

bin by depth

bin by time

decide which is best

# }
