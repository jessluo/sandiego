#
#      Analysis of the species-environment relationships through
#      regression trees
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("mvpart")
library("plyr")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())

##{ Read and reformat data ------------------------------------------------

# read data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="GMT")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards


# identify explanatory variables of interest
locVars <- c("depth", "long")
hydroVars <- c("temp", "salinity", "fluoro", "oxygen")
vars <- c(locVars, hydroVars)

# define new groups for the analysis
d$group2 <- d$group
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  return(data.frame(x[1,vars], concentration=tot))
}, .parallel=TRUE)

ggplot(dg[dg$concentration > 0,]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group2, scale="free")
ggplot(dg[dg$concentration > 0,]) + geom_histogram(aes(x=log(concentration))) + facet_wrap(~group2, scale="free")

# }


##{ Per group trees -------------------------------------------------------

m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("none"))
m <- mvpart(concentration ~ depth + long, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)

# }


##{ Community structure ---------------------------------------------------

# let's look at the influence of the different sizes of Solmaris on each other
dSol <- d[d$group=="Solmaris",]
dSol <- dcast(dSol, dateTimeB ~ taxon, value.var="concentration")
dSol <- join(dSol, unique(d[,c("dateTimeB", vars)]))

m <- mvpart(sol_small ~ temp + salinity + fluoro + oxygen + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)
m <- mvpart(sol_small ~ depth + long + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)

# }

