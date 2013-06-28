#
#      Analysis of the species-environment relationships through
#      ordination
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
  timeavg <- mean(x$dateTime)
  return(data.frame(x[1,vars], concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# casting data into wide format with group 2
dgC <- dcast (dg, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + fluoro + oxygen ~ group2, value.var="concentration")

# }


##{ Unconstrained ordination - CA  ----------------------------------------

# Correspondence Analysis
# is it going to be better to transform the concentration values? probably will be fine on either transformed values or raw values. can do a sqrt transformation to dampen the effect of dominant species. or log transformation, which will transform the concentrations (when it is not zero) to being close to normal

# start with CA with no transformation to see how it performs
dCspp <- dgC[,names(dgC) %in% c("appendicularians", "Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# CA requires all row sums to be >0 in the community data matrix
# check to see which row sums = 0 and remove them
which(rowSums(dCspp) == 0)
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
# remove the NAs
dCspp <- na.omit(dCspp)
# perform the CA using cca() from the vegan library
allCA <- cca(dCspp)
# check the summary
head(summary(allCA)) # first two CA axes explain 56% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")

# # removing the appendicularians from the CA
# dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# dCspp <- dCspp[-which(rowSums(dCspp)==0),]
# dCspp <- na.omit(dCspp)
# allCA <- cca(dCspp)
# head(summary(allCA)) # first two CA axes explain 59% of the variance
# plot(allCA, scaling=2, main="CA biplot of species concentrations")
# text(allCA, dis="sp", col="red")
# # --> removing the appendicularians improves the proportion explained of the first CA a little bit. But in this case you should not remove the appendicularians from the data because the differences that you would see in the appendicularians compared to the other data would show up in the CA biplot

# it would be nice to plot the CA plots with different colors representing different depths or longitudes. does autoplot work for CA?
# can do it by hand, just grab the relevant variables from the CA output and put it in a data frame for ggplot to plot


# CA with log transformation - this is preferred because the variables are not normally distributed. And even if you can normally distribute the abundances, with the zeros it is not normally distributed. 
# Sakina said that it would be better to perform the log transformation and perform the CA on the log transformed data anyway, even if the zeros make it so that your data is no longer normal
dCspp <- dgC[,names(dgC) %in% c("appendicularians", "doliolids", "Ctenophores", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# log transform
dCspp <- log(dCspp+1)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 65% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# the first 3 CA axes are important. Can create a ggplot just like the CA biplot and change the colors of the site points to be related to the scores along the third axis. That way you can look at the influence of the three axes at all once. 

# }

## { CA on each taxon separately ----------------------------------- 
library("pastecs")

# Ctenophores
# subsetting and casting data into wide format
ds <- d[d$group=="Ctenophores",]
dsC <- dcast (ds, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + swRho + fluoro + oxygen ~ taxon, value.var="concentration")
dCspp <- dsC[,names(dsC) %in% c("Beroida", "Bolinopsis", "Charistephane", "Dryodora glandiformis", "Haeckelia beehlri", "Hormiphora californiensis", "Juvenile Lobata", "Larval Lobata", "Mertensid", "Ocyropsis maculata", "Pleurobrachia", "Thalassocalycidae inconstans", "Velamen")]
colSums (dCspp)
# --> Pleurobrachia, Charistephane, Bolinopsis and Dryodora are the least common species

# select some taxa using Escoufier’s equivalent vectors -- escouf on dCspp 
ct.escouf <- escouf(dCspp)
plot(ct.escouf)
# --> these taxa are all different, just arbitrarily select some
# mmyLevel <- identify(ct.escouf)
# selectedTaxa <- extract(ct.escouf, level=mmyLevel)

dCspp <- dsC[,names(dsC) %in% c("Beroida", "Haeckelia beehlri", "Hormiphora californiensis", "Juvenile Lobata", "Larval Lobata", "Mertensid", "Ocyropsis maculata", "Thalassocalycidae inconstans", "Velamen")]
# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# log transform
dCspp <- log(dCspp+1)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 31% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# axes 2 and 3
plot(allCA, scaling=2, main="CA biplot of species concentrations", choice=2:3)
text(allCA, dis="sp", col="red", choice=2:3)
#

# cluster analysis on ordination results
allCA$CA$v.eig # species scores
# pick number of axes
axes <- 3
CAaxes <- allCA$CA$v.eig
CAaxes <- CAaxes[1:nrow(CAaxes),1:axes]

CAdist <- dist(CAaxes, method="euclidean")

CAclust <- hclust(CAdist, method="ward")
plot(CAclust, labels=dimnames(CAaxes)[[1]])

# Hydromedusae
# subsetting and casting data into wide format
ds <- d[d$group=="Hydromedusae",]
dsC <- dcast (ds, dateTimeB + dateTime + transect + cast + down.up + depth + long + temp + salinity + swRho + fluoro + oxygen ~ taxon, value.var="concentration")
dCspp <- dsC[,names(dsC) %in% c("Annatiara", "h1", "h10_Pegantha", "h11_Haliscera", "h13", "h15", "h2_Haliscera", "h3_Cunina", "h5_Liriope", "h5b", "h6_Solmundella", "h7_Pegantha",  "h7_Rhopalonema", "h9_Aglaura", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia", "vsh")]
colSums (dCspp)
# select some taxa
# Compute Escoufier’s equivalent vectors using escouf on dCspp 
hydro.escouf <- escouf(dCspp)
plot(hydro.escouf)
# --> again, all the taxa are different

# mmyLevel <- identify(hydro.escouf)
# selectedTaxa <- extract(hydro.escouf, level=mmyLevel)

dCspp <- dsC[,names(dsC) %in% c("Annatiara", "h1", "h10_Pegantha", "h11_Haliscera", "h15", "h2_Haliscera", "h3_Cunina", "h5_Liriope", "h5b", "h6_Solmundella", "h7_Pegantha",  "h7_Rhopalonema", "h9_Aglaura", "r3", "r4_Aegina", "r5_Eutonia", "vsh")]

# removes rows that sum to zero and are also NAs
dCspp <- dCspp[-which(rowSums(dCspp)==0),]
dCspp <- na.omit(dCspp)
# log transform
dCspp <- log1p(dCspp)
allCA <- cca(dCspp)
head(summary(allCA)) # first two CA axes explain 29% of the variance
plot(allCA, scaling=2, main="CA biplot of species concentrations")
text(allCA, dis="sp", col="red")
# axes 2 and 3
plot(allCA, scaling=2, main="CA biplot of species concentrations", choice=2:3)
text(allCA, dis="sp", col="red", choice=2:3)

# cluster analysis on ordination results
allCA$CA$v.eig # species scores
# pick number of axes
axes <- 4
CAaxes <- allCA$CA$v.eig
CAaxes <- CAaxes[1:nrow(CAaxes),1:axes]

CAdist <- dist(CAaxes, method="euclidean")

CAclust <- hclust(CAdist, method="ward")
plot(CAclust, labels=dimnames(CAaxes)[[1]])
rect.hclust(CAclust, k=4, border="red")
# --> with 4 axes you get three groups: 1) h7-pegantha, h6-solmundella, h5-liriope, h5b, h15, vsh. 2) r3, h3-cunina, h2-haliscera, h7-rhopalonema. 3) h9-aglaura, h11-haliscera, r5-eutonina, h1, h10-pegantha, and then annatiara and r4-aegina are very different
# --> with 3 CA axes you get two major groups, 1) h7-pegantha, h15, vsh, h6-solmundella, h5-liriope, h5b, h1, h10-pegantha. 2) r3, h9-aglaura, h2-haliscera, h7-rhopalonema, h3-cunina, r5-eutonia, and h11-haliscera. and annatiara and aegina are different.

# order
order <- c("Annatiara", "r4_Aegina", "h7_Pegantha", "h6_Solmundella", "h5_Liriope", "h5b", "h15", "vsh", "r3", "h3_Cunina", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura", "h11_Haliscera", "r5_Eutonia", "h1", "h10_Pegantha")

dh <- d[d$group=="Hydromedusae" & d$taxon %in% order,]
dh$taxon <- factor(dh$taxon, levels=order)

ggplot(dh[dh$concentration>0,]) + geom_violin(aes(x=taxon, y=-depth, weight=concentration, colour=taxon), alpha=0.7, scale="width") + labs(colour="Taxa", y="Depth (m)", x="Taxa")
# }

##{ Unconstrained Ordination: CA with all taxa, using results from previous ordination results ------------

# }

##{ Constrained ordination: CCA ------------------------------------------------
# which is better, RDA or CCA?
# Lots of zeros in data, CCA would be better in this case.

# setting environmental vars
# use all your environmental variables in the analysis here
dCenv <- dgC[, vars] 
# setting species variables
dCspp <- dgC[,names(dgC) %in% c("appendicularians", "Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
keep <- which(rowSums(dCspp) > 0)
dCspp <- dCspp[keep,]
dCenv <- dCenv[keep,]
# check again
nrow(dCenv) == nrow(dCspp)
dCspp <- log(dCspp + 1) # should keep log transformation here

# another option is to plot the relationship between your species abundances and environmental variables and transform the environmental variables
# if you have a polynomial relationship between species abundances and one enviornmental variable, you can code in the relationship like thus:
# dCenv$temp2 <- (dCenv$temp)^2
# in this case, it is not clear that there is a polynomial relationship (the relationship is more bimodal so you would not be able to do it in thise case.)


allCCA <- cca(dCspp ~ ., data=dCenv)
head(summary(allCCA))
# the CCA component scores are much lower than the CA component scores. and constrained proportion of inertia is 11% whereas unconstrained is 89%. this isn't that great ... but we will proceed
plot(allCCA, display=c("sp","lc","cn"), main="cca triplot conc ~ env") 
text(allCCA, dis="sp", col="red")

# test for significance of CCA results
anova(allCCA, steps=1000) # surprisingly it is significant
anova(allCCA, steps=1000, by="axis") # the first two axes are significant

# test VIFs of environmental variables - high values indicate collinearity
vif.cca(allCCA)

# compute a forward selection of explanatory variables for the CCA
ccaordistep <- ordistep(cca(dCspp ~1, data=dCenv),scope=formula(allCCA),direction="forward",pstep=1000)
# this computation indicates that you can keep all your explanatory variables

# Sakina says that you should keep all your data, if the appendicularians are very different then you would see it separated out in the CCA triplot

# }


##{ Constrained ordination: RDA ------------------------------------------------

# setting environmental vars and species vars
dCenv <- dgC[, hydroVars]
dCspp <- dgC[,names(dgC) %in% c("Ctenophores", "doliolids", "Hydromedusae", "Siphonophores", "sol_large", "sol_small")]
keep <- which(rowSums(dCspp) > 0)
dCspp <- dCspp[keep,]
dCenv <- dCenv[keep,]
# check again
nrow(dCenv) == nrow(dCspp)

# transform species concentrations
# choosing log transformation. but other ones may work as well (hellinger? what is best here?)
dCspp <- log(dCspp + 1)

# do I also need to transform env vars?
# all of the env vars are bimodal to some degree

allrda <- rda(dCspp ~ ., data=dCenv)
head(summary(allrda)) # not that bad

# adjusted R squared
RsquareAdj(allrda)
r2adj <- RsquareAdj(allrda)$adj.r.squared

# RDA triplot
plot(allrda, scaling=2, title="RDA triplot, log(concentrations) ~ env")
text(allrda, dis="sp", col="red")

# permutation tests
# testing overall significance
anova.cca(allrda, step=1000)
# testing each axis
anova.cca(allrda, step=1000, by="axis") # 4th axis not significant

# VIF scores
vif.cca(allrda) # oxygen probably needs to be removed

forward.sel(dCspp, dCenv, adjR2thresh=r2adj)

# removing oxygen because it covaries with temperature
allrda <- rda(dCspp ~ temp + fluoro + salinity, data=dCenv)
head(summary(allrda)) # not that bad

# adjusted R squared
RsquareAdj(allrda)
r2adj <- RsquareAdj(allrda)$adj.r.squared

# RDA triplot
plot(allrda, scaling=2, main="RDA triplot, log(concentrations) ~ env")
text(allrda, dis="sp", col="red")

# permutation tests
# testing overall significance
anova.cca(allrda, step=1000)
# testing each axis
anova.cca(allrda, step=1000, by="axis") # 4th axis not significant

# VIF scores
vif.cca(allrda) # oxygen probably needs to be removed

forward.sel(dCspp, dCenv, adjR2thresh=r2adj)

# }
