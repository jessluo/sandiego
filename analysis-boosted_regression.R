# boosted regression trees

library(ggplot2)
library(reshape2)
library(plyr)
library(gbm)
library(dismo)
library(gridExtra)

d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

d$front <- factor(d$front, levels=c("west", "front", "east"))

## { Boosted Regression Trees with gbm package -------------------------------
# function for R^2 for gbm
gbm.r2 <- function(gbm.obj){
	# calculate a pseudo R^2 value for a gbm model
	
	# original and fitted values
	y <- gbm.obj$data$y
	fitted <- predict.gbm(gbm.obj)
	
	# number of predictors & degrees of freedom
	n <- dim(gbm.obj$data$x.order)[1]
	p <- dim(gbm.obj$data$x.order)[2]
	df <- n - p - 1
	
	# mean squared error from fitted values
	MSE <- sum((y-fitted)^2)/df
	
	# mean sum of squares total
	MST <- var(y)
	
	# R2-adjusted for fitted values
	R2 <- 1 - (MSE / MST)
	
	return(R2)
}

gbm.results <- function(gbm.obj){
  # function call
  max.trees <- gbm.obj$n.trees
  shrinkage <- gbm.obj$shrinkage
  
  best.iter <- gbm.perf(gbm.obj, method="cv", plot.it=F)
  
  summary <- summary(gbm.obj, n.trees=best.iter, plotit=F)
  
  r2 <- gbm.r2(gbm.obj)
  
  outdf <- data.frame(max.trees, shrinkage, best.iter, r2, t(summary))[2,]
  row.names(outdf) <- 1
  
  out <- list(model_call=gbm.obj, max_trees=max.trees, shrinkage=shrinkage, best_iteration=best.iter, relative_influence=summary, r2_value=r2)
  
  return(out)
}

gbm.resultsdf <- function(gbm.obj){
  # function call
  shrinkage <- gbm.obj$shrinkage
  
  best.iter <- gbm.perf(gbm.obj, method="cv", plot.it=F)
  
  summary <- summary(gbm.obj, n.trees=best.iter, plotit=F)
  
  r2 <- gbm.r2(gbm.obj)
  
  out <- data.frame(shrinkage, best.iter, r2, t(summary))[2,]
  
  return(out)
}

?gbm
# n.trees - want something in the 1000's
# the more trees you have the more variance you can explain
# cv.folds - cross validation; 
# interaction depth - how big each tree is. the higher this number is, the more complex interactions you can fit. the higher the interaction depth and the greater the learning rate the steeper the drop off of the CV curve. BUT you want higher the number of trees, so you would want lower interaction rate and lower shrinkage. 
# do you care about autocorrelation/independence of samples? no because it is dealt w/ the CV
# bag.fraction - take the default (50%)
# train.fraction or nTrain - this is when you have a training set

##{ Hydromedusae -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="Solmaris",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.02,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Solmaris.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# # testing the dismo package
# library(dismo)
# 
# gbm_dismo <- gbm.step(d[d$taxon=="sol_large",], gbm.x=c(3,6,9:10,12:13), 
# 				gbm.y=25, 
# 				bag.fraction=0.5, 
# 				tree.complexity=3,
# 				learning.rate=0.005, 
# 				max.trees=10000, 
# 				family='gaussian')

# liriope tetraphylla
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="h5_Liriope",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.03,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Liriope.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h5_Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_vsh.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_h7Pegantha_wofront.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h7_Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_h15.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_Solmundella.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h6_Solmundella", gbm.resultsdf(m)))

# # h9_Aglaura
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="h9_Aglaura",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_h9Aglaura.txt")
# gbm.results(m)
# sink()
# 
# 
# # h11_Haliscera
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="h11_Haliscera",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.0005,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_h11Haliscera.txt")
# gbm.results(m)
# sink()
# 
# h3_Cunina
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_h3Cunina.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h3_Cunina", gbm.resultsdf(m)))

# 
# # h2_Haliscera
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="h2_Haliscera",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.0005,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_h2Haliscera.txt")
# gbm.results(m)
# sink()
# 
# # h7_Rhopalonema
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="h7_Rhopalonema",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.0005,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_h7Rhopalonema.txt")
# gbm.results(m)
# sink()

# } 

## { Apps and Doliolids --------------------------------------

# Appendicularians
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_appendicularians.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_doliolids.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# }

## { Siphonophores -------------------------------------------------
# sphaeronectes
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="Sphaeronectes",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.005,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Sphaeronectes.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="lemu",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.002,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Diphyidae(lemu).txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_MUAT.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_NABI.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_AGEL.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# # Prayidae
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Prayidae",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_Prayidae.txt")
# gbm.results(m)
# sink()
# 
# # Lilyopsis
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Lilyopsis",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_Lilyopsis.txt")
# gbm.results(m)
# sink()

# }

## { Ctenophores ----------------------------------------------------

# Ocyropsis
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="Ocyropsis maculata",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.001,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Ocyropsis.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# # thalasso
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Thalassocalycidae inconstans",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_Thalassocalycidae.txt")
# gbm.results(m)
# sink()
# 
# Larval Lobata
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Larval Lobata",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=4)
# 
# sink("plots/gbm/gbm_LarvalLobata.txt")
# gbm.results(m)
# sink()
# 
# results <- rbind(results,data.frame(taxa="Larval Lobate", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/gbm_Beroida.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# 
# # Mertensid
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Mertensid",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_Mertensid.txt")
# gbm.results(m)
# sink()
# 
# # Haeckelia beehlri
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Haeckelia beehlri",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=5)
# 
# sink("plots/gbm/gbm_Haeckelia.txt")
# gbm.results(m)
# sink()
# 
# # Velamen
# m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
#          data=d[d$taxon=="Velamen",],
#          distribution="gaussian",
#          n.trees=10000,
#          shrinkage=0.001,
#          interaction.depth=3,
#          bag.fraction=0.5,
#          cv.folds=3)
# 
# sink("plots/gbm/gbm_Velamen.txt")
# gbm.results(m)
# sink()
# 

# hormiphora
m <- gbm(concentration~front+temp+salinity+fluoro+oxygen+depth,
               data=d[d$taxon=="Hormiphora californiensis",],
               distribution="gaussian",
               n.trees=10000,
               shrinkage=0.001,
               interaction.depth=3,
               bag.fraction=0.5,
               cv.folds=5)

sink("plots/gbm/gbm_Hormiphora.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))


# }

##{ GBM Without Front ---------------------------------------------------------------

results <- data.frame()

#  solmaris
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=11000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Solmaris.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Liriope.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h5_Liriope", gbm.resultsdf(m)))

# vsh
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_vsh.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))

# h7_pegantha
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_h7Pegantha.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h7_Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_h15.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Solmundella.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h6_Solmundella", gbm.resultsdf(m)))


# h3_Cunina
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_h3Cunina.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h3_Cunina", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_appendicularians.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/regbm_doliolids.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))


# sphaeronectes
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Sphaeronectes.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Diphyidae(lemu).txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_MUAT.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_NABI.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_AGEL.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Ocyropsis.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Beroida.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~temp+salinity+fluoro+oxygen+depth,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink("plots/gbm/nofront/gbm_Hormiphora.txt")
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, "data/gbm_results.csv", row.names=F)

# }

## { Manipulate & analyze results ----------------------------------------------------
results <- read.csv("data/gbm_results.csv", stringsAsFactors=FALSE)

levels <-  c("Pegantha", "Apps", "h15", "Solmaris", "vsh", "Solmun.", "Liriope", 
           "doliolids", "Hormiphora", "Agalma", "Muggiaea", "Nanomia", "Sphaero.", 
           "Beroid", "LarvLob", "Ocyropsis", "Diphyidae", "Cunina")
labels <- c("Peg.", "Apps", "h15", "SORH", "vsh",  "SOBI", "LITE", 
            "doliolids", "HOCA", "AGEL", "MUAT", "NABI", "Sphaero.", 
            "Beroid", "LarvLob", "OCMA", "Diphy.", "Sol sp.2")


results$taxa <- factor(results$taxa, levels=levels, labels=labels)
resultsM <- melt(results, id.vars="taxa", measure.vars=c("depth", "temp", "oxygen", "fluoro", "salinity", "front"))
resultsM$value <- as.numeric(resultsM$value)

exclude <- c("LarvLob")
p <- ggplot(data=resultsM[resultsM$taxa %ni% exclude,])+ geom_bar(aes(y=value, x=variable, fill=value), stat="identity") + facet_grid(taxa~.) + scale_fill_gradient("% influence", low="grey80", high="black") + scale_x_discrete(labels=c("temp" = "T", "depth"= "D", "fluoro" = "F", "oxygen"="O", "salinity" = "S", "front"="Fr")) + scale_y_continuous(breaks=c(20, 40)) + labs(y="Relative influence (%)") + theme_bw()

pdf("plots/gbm/gbm_skinny_new.pdf", width=4.5, height=11)
print(p)
dev.off()
# }
