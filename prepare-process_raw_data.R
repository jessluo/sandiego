#
#      Read and preprocess ISIIS data
#      - reformat dates, lat/long
#      - remove erroneous data
#      - detect casts
#      ...
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("plyr")
library("stringr")
library("reshape2")
library("pastecs")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

# create the final data repository
dir.create("data", showWarnings=FALSE)

message("Read and process physical data")
##{ Read and reformat physical data ---------------------------------------

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
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="GMT"))
  # NB: we say it is GMT when it is in fact local time, just to avoid having to deal with time zones
  
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
  # we are in the western hemisphere so longitude should be negative
  d$long <- -d$long
    
  # columns that are all zero are not possible
  # they are actually missing data
  # detect them
  totCol <- colSums(d[llply(d, class) == "numeric"])
  allZeroCols <- names(totCol)[which(totCol == 0)]
  # replace the content with NA
  d[,allZeroCols] <- NA
  
  # rename some columns
  d <- rename(d, c("horizontal.vel.in.water"="horizontal.vel",
                   "irrandiance"="irradiance"
  ))

  # keep only interesting data
  d <- d[,c("transect", "dateTime", "depth", "lat", "long", "temp", "salinity", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged")]
  
  return(d)
  
}, .progress="text")

# remove adply crap
phy <- phy[,-1]

# }


##{ Check and correct the physical data -----------------------------------

summary(phy)

# inspect water mass data
phyM <- melt(phy, id.vars=c("dateTime"), measure.vars=c("depth", "temp", "salinity", "fluoro", "oxygen", "irradiance"))
ggplot(data=phyM) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scales="free")

# remove erroneous points when ISIIS lost contact, in transect 2
# NB: this is specific to this dataset, a general solution would be better
ggplot(data=phy) + geom_point(aes(x=long, y=-depth), size=0.5) + facet_grid(transect ~.) 
which(abs(diff(phy$depth)) > 20)
phy[49417:49428,]
phy[49418:49427,]
phy <- phy[-c(49418:49427),]

# the depth gets stuck from time to time and that results in jumps afterwards. Remove those stuck points and reinterpolate the depth linarly using time.
# assign the depths in which the difference before the previous depth is 0 to be NA
phy$depth[which(diff(phy$depth)==0)+1] <- NA

# interpolation
phy$depth <- approx(phy$dateTime, phy$depth, phy$dateTime, method="linear")$y

# irradiance looks funny, all negative
ggplot(data=phy) + geom_histogram(aes(x=irradiance)) + facet_wrap(~transect)
# same for all transects

# look at profiles
ggplot(phy) + geom_path(aes(x=temp, y=-depth), alpha=0.5) + facet_wrap(~transect)
ggplot(phy) + geom_path(aes(x=salinity, y=-depth), alpha=0.5) + facet_wrap(~transect)
ggplot(phy) + geom_path(aes(x=fluoro, y=-depth), alpha=0.5) + facet_wrap(~transect)
ggplot(phy) + geom_path(aes(x=oxygen, y=-depth), alpha=0.5) + facet_wrap(~transect)
ggplot(phy) + geom_path(aes(x=irradiance, y=-depth), alpha=0.5) + facet_wrap(~transect)
ggplot(phy) + geom_path(aes(x=irradiance, y=-depth), alpha=0.5) + facet_wrap(~transect) + scale_x_continuous(limits=c(-1.7E-6, -7.5E-7))
# real issue with irradiance
# the sensor was actually broken, remove all values
phy$irradiance <- NA


# inspect navigation data
ggplot(data=phy) + geom_histogram(aes(x=horizontal.vel), binwidth=100)
ggplot(data=phy) + geom_histogram(aes(x=abs(vertical.vel)), binwidth=100)
ggplot(data=phy) + geom_histogram(aes(x=abs(pitch)))

# obvious issues with horizontal velocities
length(which(phy$horizontal.vel <= 0))
ggplot(data=phy[which(phy$horizontal.vel <= 2000),]) + geom_point(aes(x=long, y=-depth, size=-horizontal.vel, colour=horizontal.vel <= 0)) + facet_grid(transect ~.)
# ggplot(data=phy) + geom_point(aes(x=long, y=-depth, colour=horizontal.vel, size=abs(horizontal.vel))) + facet_grid(transect ~.) + scale_colour_gradient2()
ggplot(phy) + geom_path(aes(x=horizontal.vel, y=-depth), alpha=0.5) + facet_wrap(~transect)
# remove everything around the main peak
phy$horizontal.vel[phy$horizontal.vel <= 1500 | phy$horizontal.vel >= 3200] <- NA

# issues with vertical vel too
ggplot(phy) + geom_path(aes(x=vertical.vel, y=-depth), alpha=0.5) + facet_wrap(~transect)
# look at it in time
phyT2 <- phy[phy$transect==2,]
ggplot(phyT2[2500:5000,]) + geom_path(aes(x=dateTime, y=vertical.vel, colour=depth))
# compute percentage of "extreme" values
sum(abs(phy$vertical.vel) > 1000, na.rm=T) / length(na.omit(phy$vertical.vel))
# [1] 0.02474908
# that's actually noise because of a power outage for the ADCP

# recompute the vertical velocity from the depth change / time
phy$vertical.vel <- c(NA, diff(phy$depth) / as.numeric(diff(phy$dateTime))) * 1000
ggplot(data=phy) + geom_histogram(aes(x=vertical.vel), binwidth=100)
ggplot(phy) + geom_path(aes(x=vertical.vel, y=-depth), alpha=0.5) + facet_wrap(~transect)
# still some peaks, far less

# remove peaks
phy$vertical.vel[abs(phy$vertical.vel) > 1000] <- NA


# calculate in water velocity
phy$velocity <- sqrt(phy$horizontal.vel^2 + phy$vertical.vel^2)
# convert it in m/s
phy$velocity <- phy$velocity / 1000

# inspect
ggplot(data=phy) + geom_histogram(aes(x=velocity), binwidth=0.01)
# mean velocity is really close to 2.5, which is what we want
mean(phy$velocity, na.rm=TRUE)

# }


##{ Detect up and down casts and number them ------------------------------

# there is still an issue here with the numbering of the casts. in transect 2, the first minute and a half of data is at the surface and it is marked as an upcast. 
# brute force method of removal -- as a temporary fix so I can do the joining. think of a better way to detect this and remove the data
phy <- phy[-which(phy$transect==2)[1:400],]

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


# check cast numbering
ggplot(data=phy) + geom_path(aes(x=long, y=-depth, colour=as.factor(cast), linetype=factor(down.up))) + facet_grid(transect ~.) 

# }

# save it as text
write.csv(phy, "data/phy.csv", row.names=FALSE)



message("Read and process biological data")

##{ Biological data -------------------------------------------------------

bioData <- "raw_biological_data"

# NB: The data files for various groups have to be read separately because the format and convention change between each group

# SOLMARIS
solFiles <- list.files(bioData, pattern="Solmaris", full=TRUE)
sol <- adply(solFiles, 1, read.csv, stringsAsFactors=FALSE)

# identify transects (based on the day number)
sol$X1 <- solFiles[sol$X1]
sol$X1 <- str_replace(sol$X1, bioData, "")
sol$transect <- as.numeric(str_sub(sol$X1, 4, 5)) - 14

# create a true date+time column
sol$dateTime <- as.POSIXct(str_c(sol$date, " ", sprintf("%02i",sol$hour), ":", sprintf("%02i",sol$min), ":", sprintf("%02i",sol$sec), ".", sol$s.1000), tz="GMT")

# convert to the tall format
sol <- sol[ , ! names(sol) %in% c("X1", "date", "hour", "min", "sec", "s.1000")]
solT <- melt(sol, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
solT$taxon <- as.character(solT$taxon)
solT$group <- "Solmaris"
# remove the zeros
solT <- solT[solT$count != 0,]

# add subsampling interval (Solmaris were looked for in one every five frame)
solT$sub <- 5


# SIPHONOPHORES
siphFiles <- list.files(bioData, pattern="Siphonophores", full=TRUE)
siph <- adply(siphFiles, 1, read.csv, stringsAsFactors=FALSE)

# identify transects (based on the day number)
siph$X1 <- siphFiles[siph$X1]
siph$X1 <- str_replace(siph$X1, bioData, "")
siph$transect <- as.numeric(str_sub(siph$X1, 4, 5)) - 14

# create a true date+time column
siph$dateTime <- as.POSIXct(str_c(siph$date, " ", sprintf("%02i",siph$hour), ":", sprintf("%02i",siph$min), ":", sprintf("%02i",siph$sec), ".", siph$s.1000), tz="GMT")

# combine with and without tail
siph$Type1 <- siph$Type1 + siph$Type1_wotail
siph$Type2 <- siph$Type2 + siph$Type2_wotail
siph$Type3 <- siph$Type3 + siph$Type3_wotail
siph <- siph[,-which(str_detect(names(siph), "_wotail"))]

# give siphonophore taxa meaningful names
siph <- rename(siph, c("Type1"="Diphyidae", "Type2"="Sphaeronectes", "Type3"="Physonect", "Type4"="Prayidae", "Type5"="Lilyopsis"))

# convert to the tall format
siph <- siph[ , ! names(siph) %in% c("X1", "date", "hour", "min", "sec", "s.1000")]
siphT <- melt(siph, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
siphT$taxon <- as.character(siphT$taxon)
siphT$group <- "Siphonophores"
# remove the zeros
siphT <- siphT[siphT$count != 0,]

# add subsampling interval (Siphonophores were looked for in every frame)
siphT$sub <- 1


# CTENOPHORES
cteFiles <- list.files(bioData, pattern="Ctenophores", full=TRUE)
cteT <- adply(cteFiles, 1, function(file){
  d <- read.csv(file, stringsAsFactors=FALSE)
  
  # clean up taxon names
  d$species <- str_trim(d$species)
  
  # recompute the time
  # the record contains the time of the first frame in the stack and the frame number within the stack
  # there are 300 frames in a stack and they represent a time of 17.375 seconds
  d$dateTime <- as.POSIXct(str_c(d$year, "-", sprintf("%02i",d$month), "-", sprintf("%02i",d$date), " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec)), tz="GMT")
  d$dateTime <- d$dateTime + 17.375 * d$x.300 / 300
  
  # extract downcast number
  downcast <- d$folder
  downcast <- str_sub(downcast, start=4, end=5)
  downcast <- str_trim(downcast)
  downcast <- as.numeric(downcast)
  d$downcast <- downcast
  
  # add a reorder columns
  d$group <- "Ctenophores"
  d <- rename(d, c("species"="taxon"))
  d <- d[, c("downcast", "dateTime", "group", "taxon")]
  d$count <- 1
  
  return(d)
})

# correct spelling on one taxon
cteT$taxon[which(cteT$taxon=="Ocyropsis maciluta")] <- "Ocyropsis maculata"

# identify transects (based on the day number)
cteT$X1 <- cteFiles[cteT$X1]
cteT$X1 <- str_replace(cteT$X1, bioData, "")
cteT$transect <- as.numeric(str_sub(cteT$X1, 4, 5)) - 14
cteT <- cteT[,-1]

# add subsampling interval (Ctenophores were looked for in every frame)
cteT$sub <- 1


# HYDRO + DOLIOLIDS + APPS
hFiles <- list.files(bioData, pattern=glob2rx("*Hydromedusae*csv"), full=TRUE)
h <- adply(hFiles, 1, function(file) {
  d <- read.csv(file, stringsAsFactors=FALSE)

  # homogenise taxa names
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
  d$dateTime <- as.POSIXct(str_c("2010-10-", d$day, " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec), ".", d$s.1000), tz="GMT")
  # remove extraneous date columns
  d <- d[,!names(d) %in% c("day", "hour", "min", "sec", "s.1000")]
  
  return(d)
})

# remove taxa we are not interested in
h <- h[,!names(h) %in% c("muat", "lemu", "spko", "si.2.eu", "nabi", "agel", "coor")]
# and the stuff that's added by adply
h <- h[,!names(h) %in% c("X1")]

# give hydromedusae meaningful names
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

# add the group identification
hT$group <- ifelse(hT$taxon == "doliolids", "Tunicates", "Hydromedusae")

# add subsampling interval
# Hydromedusae were looked for in every frame
hT$sub <- 1
# Appendicularians were looked for in 1 fifth of a frame, every 20 frames
appT$sub <- 100
# except Appendicularians in downcast 1 of transect 1, looked for in one fifth of a frame, every 10 frames
appT$sub[which(appT$transect == 1 & appT$downcast == 1)] <- 50


# put all the data together
bio <- rbind.fill(solT, siphT, cteT, appT, hT)

# shift the time by 3 hours
bio$dateTime <- bio$dateTime - 3 * 3600

# change the name of casts and do not mix it with the cast number defined by the physics just yet
# they should be the same but we never know
bio <- rename(bio, c("downcast"="cast.bio"))

# }


##{ Check biological data -------------------------------------------------

head(bio)
unique(bio$transect)
unique(bio$cast.bio)
unique(bio$group)
unique(bio$taxon)
range(bio$dateTime)
unique(bio$count)
count(bio$count)

countPerTransect <- ddply(bio, ~transect + group + taxon, function(x) {sum(x$count, na.rm=T)})
dcast(countPerTransect, group+taxon~transect, value.var="V1")

# -> all that seems consistent and OK
# note that upuon first inspection, it seems strange that there are so much less appendicularians in transect 2

# }

# save it as text
write.csv(bio, "data/bio.csv", row.names=FALSE)
