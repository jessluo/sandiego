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

# This section performs the following tasks:
# 1) read in the physical data files from a folder
# 2) cleans up the names by removing the extra characters, units and making everything lowercase
# 3) from the time format in the physical data, converts it to the universal time format
# 4) corrects the time zone and adds in transect number, and converts lat/long into decimal degrees

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
  d$dateTime <- as.POSIXct(strptime(d$dateTime, format="%m/%d/%y %H:%M:%OS", tz="America/Los_Angeles"))
  # NB: we say it is America/Los_Angeles when it is in fact local time, just to avoid having to deal with time zones
  
  # shift all physical data back 3 hours
  # NB: this is specific for the San Diego dataset, as the physical data was recorded in Eastern Time
  # and our sampling was done in Pacific Time
  d$dateTime <- d$dateTime - 3 * 3600
  
  # code in a transect number
  # this is not robust for all physical data but is necessary here
  # subtract 14 because our first sampling day was Oct 15, which is transect 1
  d$transect <- dd-14
  
  # reformat the lat and long in decimal degrees
  to.dec <- function(x) {
    # split in degree, minute, second
    pieces <- str_split_fixed(x, "° |'", 3)
    # extract orientation (S/N and E/W)
    orientation <- str_sub(pieces[,3], -1)
    # remove orientation to only keep numbers
    pieces[,3] <- str_replace(pieces[,3], "[NSEW]", "")
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
  d <- d[,c("transect", "dateTime", "depth", "lat", "long", "temp", "salinity", "pressure", "fluoro", "oxygen", "irradiance", "heading", "horizontal.vel", "vertical.vel", "pitch", "vol.imaged")]
  
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
# the horizontal velocity should be the velocity parallel to the direction of movement (the velocity of the instrument)
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


# in water velocity is the horizontal velocity
phy$velocity <- phy$horizontal.vel
# convert it in m/s
phy$velocity <- phy$velocity / 1000

# inspect
ggplot(data=phy) + geom_histogram(aes(x=velocity), binwidth=0.01)
# mean velocity is really close to 2.5, which is what we want
mean(phy$velocity, na.rm=TRUE)

phy <- phy[, -which(names(phy)=="horizontal.vel")]
# }

##{ Calculate seawater density --------------------------------------------

library(oce)

# using salinity, temperature and pressure for calculation of seawater density
# using UNESCO formulation
phy$swRho <- swRho(phy$salinity, phy$temp, phy$pressure, eos="unesco")

# not keeping pressure anymore
phy <- phy[, -which(names(phy)=="pressure")]

# }

##{ Detect up and down casts and number them ------------------------------

# there is still an issue here with the numbering of the casts. in transect 2, the first minute and a half of data is at the surface and it is marked as an upcast. 
# brute force method of removal -- as a temporary fix so I can do the joining. think of a better way to detect this and remove the data
phy$flag <- FALSE
phy$flag[which(phy$transect==2)[1:370]] <- TRUE
phy$flag[which(phy$transect==3)[1:50]] <- TRUE
# p <- ggplot(data=phy) + geom_path(aes(x=long, y=-depth, colour=flag), size=0.5) + facet_grid(transect ~.)
# p + xlim(-119.15, -119.1)
# p + xlim(-119.95, -119.9)
phy <- phy[!phy$flag,]

phy <- ddply(phy, ~transect, function(d) {
  # detect up and down casts by smoothing the depth profile and finding the turning points
  # smooth depths
  order <- 25
  # tested options between 10-20, correct number after order=14. use order between 16-20 for robustness
  depth_avg <- decaverage(-d$depth, times=3, weights=c(seq(1, order), order+1, seq(order, 1, -1)))
  # plot(depth_avg)
  depth_avg <- as.numeric(pastecs::extract(depth_avg, component="filtered"))
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
ddply(phy, ~transect+cast, nrow)

ggplot(data=phy) + geom_path(aes(x=long, y=-depth, colour=factor(cast), linetype=down.up), size=0.5) + facet_grid(transect ~.) 

# inspect turniing points to check the accuracy of the detection of up and down casts
# p <- ggplot(data=phy) + geom_point(aes(x=long, y=-depth, colour=down.up), size=1) + facet_grid(transect ~.)
# pieces <- seq(min(phy$long), max(phy$long), length.out=8)
# for (i in 1:(length(pieces)-1)) {
#   print(p + ylim(-7,0) + xlim(pieces[i], pieces[i+1]))
#   print(p + ylim(-138,-123) + xlim(pieces[i], pieces[i+1]))
# }
# it is better at the surface than at depth, but good enough

# remove the flag column
phy <- phy[, which(names(phy) != "flag")]

# delineate the front
# initialize
phy$front <- NA

# delineate the frontal region
phy[phy$transect==1 & phy$cast <=11,]$front <- "east"
phy[phy$transect==1 & phy$cast >=12 & phy$cast <= 15,]$front <- "front"
phy[phy$transect==1 & phy$cast >=16,]$front <- "west"
phy[phy$transect==2 & phy$cast <=22,]$front <- "east"
phy[phy$transect==2 & phy$cast >=23 & phy$cast <= 28,]$front <- "front"
phy[phy$transect==2 & phy$cast >=29,]$front <- "west"
phy[phy$transect==3 & phy$cast <=8,]$front <- "west"
phy[phy$transect==3 & phy$cast >=9 & phy$cast <= 13,]$front <- "front"
phy[phy$transect==3 & phy$cast >=14,]$front <- "east"
phy$front <- factor(phy$front, levels=c("west", "front", "east"))


# }


# save it as text
write.csv(phy, "data/phy.csv", row.names=FALSE)



message("Read and process biological data")

##{ Biological data -------------------------------------------------------

bioData <- "raw_biological_data"

# NB: The data files for various groups have to be read separately because the format and convention change between each group

# 
# SOLMARIS
solFiles <- list.files(bioData, pattern="Solmaris", full=TRUE)
sol <- adply(solFiles, 1, read.csv, stringsAsFactors=FALSE)

# identify transects (based on the day number)
sol$X1 <- solFiles[sol$X1]
sol$X1 <- str_replace(sol$X1, bioData, "")
sol$transect <- as.numeric(str_sub(sol$X1, 4, 5)) - 14

# create a true date+time column
sol$dateTime <- as.POSIXct(str_c(sol$date, " ", sprintf("%02i",sol$hour), ":", sprintf("%02i",sol$min), ":", sprintf("%02i",sol$sec), ".", sol$s.1000), tz="America/Los_Angeles")

# convert to the tall format
sol <- sol[ , ! names(sol) %in% c("X1", "date", "hour", "min", "sec", "s.1000")]
solT <- melt(sol, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
solT$taxon <- "Solmaris"
solT$group <- "Hydromedusae"
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
siph$day <- as.numeric(str_sub(siph$X1, 4, 5))

# shift by one day when we cross midnight
siph$day <- ifelse(siph$hour >= 18 & siph$hour <= 23, siph$day, siph$day+1)

# convert that into POSIXct
siph$dateTime <- as.POSIXct(str_c("2010-10-", siph$day, " ", sprintf("%02i",siph$hour), ":", sprintf("%02i",siph$min), ":", sprintf("%02i",siph$sec), ".", siph$s.1000), tz="America/Los_Angeles")

# combine with and without tail
siph$Type1 <- siph$Type1 + siph$Type1_wotail
siph$Type2 <- siph$Type2 + siph$Type2_wotail
siph$Type3 <- siph$Type3 + siph$Type3_wotail
siph <- siph[,-which(str_detect(names(siph), "_wotail"))]

# give siphonophore taxa meaningful names
siph <- rename(siph, c("Type1"="Diphyidae", "Type2"="Sphaeronectes", "Type3"="Physonect", "Type4"="Prayidae", "Type5"="Lilyopsis", "LEMU"="lemu", "MUAT"="muat", "NABI"="nabi", "AGEL"="agel", "COOR"="coor", "FOED"="foed"))

# convert to the tall format
siph <- siph[ , ! names(siph) %in% c("X1", "day", "hour", "min", "sec", "s.1000")]
siphT <- melt(siph, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
siphT$taxon <- as.character(siphT$taxon)
siphT$group <- "Siphonophores"

# remove the zeros and NAs
siphT <- siphT[complete.cases(siphT),]
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
  # there are 300 frames in a stack and they represent a time of 17.55 seconds
  d$dateTime <- as.POSIXct(str_c(d$year, "-", sprintf("%02i",d$month), "-", sprintf("%02i",d$date), " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec)), tz="America/Los_Angeles")
  d$dateTime <- d$dateTime + 17.55 * d$x.300 / 300
  
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
  
  # TODO: consider whether you want to split up the H7 groups in the future, as Rhopalonema and Cucina are fairly different
  # in third transect, we took note of the differences between the two and recorded them as different
  # for the purposes of having a homogonized data set, we will combine them here
  if ("h7.cucina" %in% names(d)){
    d$h7 <- d$h7 + d$h7.cucina
    d <- d[,!names(d) %in% "h7.cucina"]
  }
  
  # combine the Liriope tetraphylla
  d$h5a <- d$h5b + d$h5a
  d <- d[,!names(d) %in% "h5b"]
  
  d$h14 <- ifelse (("h14" %in% names(d)), d$h14, 0)
  
  # compute the date + time
  # detect starting day from the file name
  file <- str_replace(file, bioData, "")
  day <- as.numeric(str_sub(file, 4, 5))
  # use that to determine the transect number
  d$transect <- day - 14
  # shift by one day when we cross midnight
  d$day <- ifelse(d$hour >= 18 & d$hour <= 23, day, day+1)
  # convert that into POSIXct
  d$dateTime <- as.POSIXct(str_c("2010-10-", d$day, " ", sprintf("%02i",d$hour), ":", sprintf("%02i",d$min), ":", sprintf("%02i",d$sec), ".", d$s.1000), tz="America/Los_Angeles")
  # remove extraneous date columns
  d <- d[,!names(d) %in% c("day", "hour", "min", "sec", "s.1000")]
  
  return(d)
})

# pull out the other siphonophores
# subset the hydromedusae, apps and siphs dataframe for just the siphs
siph2 <- h[,c("downcast", "transect", "dateTime", "muat", "lemu", "nabi", "agel", "coor", "foed", "spko", "si.2.eu")]

# probably not necessary to include spko and si.2.eu in this case
siph2 <- siph2[, !names(siph2) %in% c("spko", "si.2.eu")]

# replace the rest of the NAs with zeros
siph2$lemu[is.na(siph2$lemu)] <- 0
siph2$muat[is.na(siph2$muat)] <- 0
siph2$nabi[is.na(siph2$nabi)] <- 0
siph2$foed[is.na(siph2$foed)] <- 0
siph2$agel[is.na(siph2$agel)] <- 0
siph2$coor[is.na(siph2$coor)] <- 0

# convert to tall format
siph2T <- melt(siph2, id.vars=c("transect", "downcast", "dateTime"), variable.name="taxon", value.name="count")
siph2T$taxon <- as.character(siph2T$taxon)
# remove the zeros
siph2T <- siph2T[siph2T$count != 0,]

# add subsampling rate
siph2T$sub <- 1

# add group name
siph2T$group <- "Siphonophores"


# remove taxa we are not interested in
h <- h[,!names(h) %in% c("h8", "muat", "lemu", "spko", "si.2.eu", "nabi", "agel", "coor", "foed")]
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
bio <- rbind.fill(solT, siphT, siph2T, cteT, appT, hT)
# bio <- rbind.fill(solT, siphT, cteT, appT, hT)

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
