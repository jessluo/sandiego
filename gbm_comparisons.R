

## { Set Run Permutation --------------------------------------
N = 161
dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+depth+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 162
dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~depth+front+oxygen+fluoro+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 163

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~depth+front+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 164

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+depth+fluoro+oxygen+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 165

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+depth+oxygen+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 166

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+oxygen+depth+salinity+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 167

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+oxygen+salinity+depth+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 168

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+oxygen+salinity+temp+depth,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 169

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+salinity+oxygen+temp+depth,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }


## { Set Run Permutation --------------------------------------
N = 170

dir.create(str_c("plots/gbm/COMPARISONS/", N, sep=""), showWarnings=F)

print(str_c("OLD: concentration~",str_c(perdf[N-1,1], perdf[N-1,2], perdf[N-1,3], perdf[N-1,4], perdf[N-1,5], perdf[N-1,6], sep="+"), sep=""), quote=F)

print(str_c("NEW: concentration~",str_c(perdf[N,1], perdf[N,2], perdf[N,3], perdf[N,4], perdf[N,5], perdf[N,6], sep="+"), sep=""), quote=F)


# }

##{ Run all GBMs -----------------------------------------------
results <- data.frame()

#  solmaris
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="Solmaris",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.02,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmaris.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmaris", gbm.resultsdf(m)))

# liriope tetraphylla
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="h5_Liriope",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.03,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Liriope.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Liriope", gbm.resultsdf(m)))

# vsh

m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="vsh",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.015,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_vsh.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="vsh", gbm.resultsdf(m)))


# h7_pegantha
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="h7_Pegantha",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Pegantha.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Pegantha", gbm.resultsdf(m)))


# h15
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="h15",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.009,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_h15.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="h15", gbm.resultsdf(m)))


# h6_Solmundella
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="h6_Solmundella",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Solmundella.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Solmundella", gbm.resultsdf(m)))

# h3_Cunina
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="h3_Cunina",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sol_sp2.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sol_sp2", gbm.resultsdf(m)))

# Appendicularians
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="appendicularians",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_appendicularians.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="appendicularians", gbm.resultsdf(m)))

# Doliolids
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="doliolids",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.0075,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_doliolids.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="doliolids", gbm.resultsdf(m)))

# sphaeronectes
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="Sphaeronectes",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Sphaeronectes.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Sphaeronectes", gbm.resultsdf(m)))

# diphyidae
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="lemu",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.002,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Diphyidae.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Diphyidae", gbm.resultsdf(m)))

# muggiaea
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="muat",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Muggiaea.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Muggiaea", gbm.resultsdf(m)))

# Nanomia bijuga
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="nabi",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.005,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Nanomia.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Nanomia", gbm.resultsdf(m)))

# Agalma elegans
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="agel",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.003,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Agalma.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Agalma", gbm.resultsdf(m)))

# Ocyropsis
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="Ocyropsis maculata",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Ocyropsis.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Ocyropsis", gbm.resultsdf(m)))

# Beroida
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="Beroida",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Beroida.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Beroida", gbm.resultsdf(m)))

# hormiphora
m <- gbm(concentration~front+fluoro+salinity+oxygen+depth+temp,
         data=d[d$taxon=="Hormiphora californiensis",],
         distribution="gaussian",
         n.trees=10000,
         shrinkage=0.001,
         interaction.depth=3,
         bag.fraction=0.5,
         cv.folds=5)

sink(str_c("plots/gbm/COMPARISONS/", N, "/gbm_Hormiphora.txt", sep=""))
gbm.results(m)
sink()

results <- rbind(results,data.frame(taxa="Hormiphora", gbm.resultsdf(m)))

# }

##{ Print and save results ----------------------------------------------------------
row.names(results) <- 1:nrow(results)

write.csv(results, str_c("plots/gbm/COMPARISONS/gbm_results_", N, ".csv", sep=""), row.names=F)

# }