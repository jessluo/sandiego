#
#      Analysis of the species-environment relationships through
#      regression trees
#
#  (c) Copyright 2013 Jean-Olivier Irisson
#                     Jessica Luo
#      GNU General Public License v3
#
#--------------------------------------------------------------------------

library("vegan")
library("mvpart")
library("packfor")
library("plyr")
library("reshape2")
library("ggplot2")
library("foreach")
library("doParallel")
registerDoParallel(cores=detectCores())

##{ Read and reformat data ------------------------------------------------

# read data
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")
# NB: make sure time is set in GMT (even if it wasn't) to avoid dealing with tz afterwards

# rename env var names
d <- rename(d, replace=c("salinity" = "sal", "oxygen" = "oxy", "fluoro" = "flu", "swRho" = "rho"))

# identify explanatory variables of interest
locVars <- c("depth", "long", "front")
hydroVars <- c("temp", "sal", "flu", "oxy")
vars <- c(locVars, hydroVars)

# define broad groups for the analysis
d$group2 <- d$group
d$group2[d$group == "Solmaris"] <- d$taxon[d$group == "Solmaris"]
d$group2[d$group == "Tunicates"] <- d$taxon[d$group == "Tunicates"]
d[d$group2=="appendicularians","group2"] <- "Appendicularians"
d[d$group2=="doliolids","group2"] <- "Doliolids"

# compute total concentration per group
dg <- ddply(d, ~transect + cast + down.up + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# visualize
# ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=concentration)) + facet_wrap(~group2, scale="free")
# ggplot(dg[which(dg$concentration > 0),]) + geom_histogram(aes(x=log(concentration))) + facet_wrap(~group2, scale="free")

# casing data into wide format with group 2
dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + sal + flu + oxy ~ group2, value.var="concentration")

# }

## { Define secondary groups for analysis -------------------------------
# use the ordination results to define groups for analysis
# hydromedusae: excluding the species that are present in less than 1% of the population (and keeping broad taxonomic groups), these are the following groups:
# 1) deep trachymedusae: h2-haliscera, h7-rhopalonema, h9-aglaura, h11-haliscera, 2) liriope and h5b, 3) shallow narcos: h7-pegantha, h6-solmundella, 4) deep narco's -- h3-cunina, and 5) other hydros: h1, h15, vsh
# ctenophores: keep broad taxonomic groups: cydippids, lobates, beroe, thalasso, velamen
# siphs: keep groups: physonectae, diphyidae, sphaeronectes, prayidae

# select downcasts for inclusion
dcchoose <- unique(d[d$group=="Hydromedusae" & d$transect==1,"cast"])

# define functional/taxonomic groups
Deep_Trachy <- c("h11_Haliscera", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura")
Shallow_Trachy <- c("h5_Liriope", "h5b")
Shallow_Narco <- c("h7_Pegantha", "h6_Solmundella")
Deep_Narco <- c("h3_Cunina")
Other_Hydro <- c("h1", "h15", "vsh")
Cydippida <- c("Haeckelia beehlri", "Hormiphora californiensis", "Mertensid")
Lobata <- c("Bolinopsis", "Ocyropsis maculata", "Juvenile Lobata", "Larval Lobata")
Prayidae <- c("Lilyopsis", "Prayidae")

# assign these groups into a different column
d$group3 <- d$taxon
d$group3[d$taxon %in% Deep_Trachy] <- "Deep Trachy"
d$group3[d$taxon %in% Shallow_Trachy] <- "Shallow Trachy"
d$group3[d$taxon %in% Shallow_Narco] <- "Shallow Narco"
d$group3[d$taxon %in% Deep_Narco] <- "Deep Narco"
d$group3[d$taxon %in% Other_Hydro] <- "Other Hydro"
d$group3[d$taxon %in% Cydippida] <- "Cydippida"
d$group3[d$taxon %in% Lobata] <- "Lobata"
d$group3[d$taxon %in% Prayidae] <- "Prayidae"

# exclude rare taxa
`%ni%` <- Negate(`%in%`) 

exclude <- c("Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown", "Annatiara", "h10_Pegantha", "h13", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia")
d <- d[d$group3 %ni% exclude,]

# select which portions of the water column to sample
drt <- d[d$transect==1 & d$cast %in% dcchoose,]
drt <- rbind(drt, d[d$transect==2,])

# compute total concentration per group
drt <- ddply(drt, ~transect + cast + front + dateTimeB + group3, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# }

##{ Creating some species plots ------------------------------------------
dO=d[d$taxon=="Ocyropsis maculata",]
dO$lconc=log1p(dO$concentration)
dO$pres <- dO$concentration > 0

ggplot(dO) + geom_point(aes(x=salinity, y=temp, size=pres, shape=long>119.6, colour=depth> 50), alpha=0.7) + scale_size_manual(values=c(1, 2.5))

dO=d[d$group=="Solmaris",]
dO$lconc=log1p(dO$concentration)
dO$pres <- dO$concentration > 0

ggplot(dO) + geom_point(aes(x=salinity, y=temp, size=concentration, color=taxon, alpha=pres)) 
#}

##{ Per group trees -------------------------------------------------------

m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("none"))
m <- mvpart(concentration ~ depth + long, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="sol_small",], xv=("min"), xval=5, xvmult=100)

# apps are not very interesting
m <- mvpart(concentration ~ depth + long + temp + salinity + fluoro + oxygen, data=dg[dg$group2=="appendicularians",], xv=("min"), xval=5, xvmult=100)

# hydromedusae
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Hydromedusae",], xv=("1se"), xval=5, xvmult=100)

# siphonophores
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Siphonophores",], xv=("min"), xval=5, xvmult=100)

# ctenophores - cannot explain ctenophores with this method because they are not numerous enough. can only use logistic regression
m <- mvpart(concentration ~ temp + salinity + fluoro + oxygen, data=dg[dg$group2=="Ctenophores",], xv=("min"), xval=5, xvmult=100)


# }

##{ Community structure ---------------------------------------------------

# let's look at the influence of the different sizes of Solmaris on each other
dSol <- d[d$group=="Solmaris",]
dSol <- dcast(dSol, dateTimeB ~ taxon, value.var="concentration")
dSol <- join(dSol, unique(d[,c("dateTimeB", vars)]))

m <- mvpart(sol_small ~ temp + salinity + fluoro + oxygen + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)
m <- mvpart(sol_small ~ depth + long + sol_large, data=dSol, xv=("1se"), xval=5, xvmult=100)

# in De'ath (2002) he says to standardize the species and then sites to the same mean. And also to standardize the data to give chi-squared distances between sites. this then gives an opportunity to compare with cca
# decostand
# setting environmental vars
dCenv <- dgC[, hydroVars]
spp <- c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")
dCspp <- dgC[,names(dgC) %in% spp]
keep <- which(rowSums(dCspp) > 0)
dCspp <- dCspp[keep,]
dCenv <- dCenv[keep,]
# check again
nrow(dCenv) == nrow(dCspp)

dCbind <- cbind(dCspp, dCenv)

m <- mvpart(data.matrix(decostand(dCbind[,names(dCbind) %in% spp], method="chi.square")) ~  temp + salinity + fluoro + oxygen, data=dCbind, xv=("min"), xval=5, xvmult=100)


# }

##{ Logistic Regression  ------------------------------------------------
# use logistic regression on presence / absence data

# code in Presence or Absence in original data frame
d$pa <- 0
d$pa[which(d$concentration != 0)] <- 1


# cast data frame
dC <- dcast (d, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + fluoro + oxygen ~ taxon, value.var="pa")
# make names you can call
names(dC) <- make.names(names(dC))

# try logistic regression on Ocyropsis maculata (Ctenophore)
nrow(dC)

ocmalogit <- glm(Ocyropsis.maculata ~ (temp + salinity + fluoro + oxygen)^4 + I(temp^2) + I(salinity^2) + I(fluoro^2) + I(oxygen^2) , data=dC, family=binomial("logit"))
summary(ocmalogit) # intercept, salinity, and oxygen are significant explanatory variables
# make a model only out of significant explanatory variables

m <- update(ocmalogit, . ~ . - temp:salinity:fluoro:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:fluoro)
summary(m)
m <- update(m, . ~ . - salinity:fluoro:oxygen)
summary(m)

m <- step(ocmalogit, direction="backward")
summary(m)
m <- update(m, . ~ . - salinity:fluoro:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - salinity:fluoro)
summary(m)
m <- update(m, . ~ . - I(fluoro^2))
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)
m <- update(m, . ~ . - temp:salinity:oxygen)
summary(m)




ocmalogit <- glm(Ocyropsis.maculata ~ salinity + oxygen, data=dC, family=binomial("logit"))
summary(ocmalogit)
plot(predict(ocmalogit, type="response") ~ salinity + oxygen, data=dC)

# TODO: try others

# }

