#
#       Analysis of community structure via correlation
#
#       Jessica Luo, 2013
#
#--------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read data -------------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")

# }

## { Construct a correlation heatmap for all taxa w/ pre-defined groups ----------------------------------
Deep_Hydro <- c("h11_Haliscera", "h2_Haliscera", "h7_Rhopalonema", "h9_Aglaura", "h3_Cunina")
Liriope <- c("h5_Liriope", "h5b")
Shallow_Narco <- c("h7_Pegantha", "h6_Solmundella")
Other_Hydro <- c("h1", "h15", "vsh")
Cydippida <- c("Haeckelia beehlri", "Hormiphora californiensis", "Mertensid")
Lobata_Thalasso <- c("Ocyropsis maculata", "Juvenile Lobata", "Larval Lobata", "Thalassocalycidae inconstans")
Prayidae <- c("Lilyopsis", "Prayidae")
Solmaris <- c("sol_large", "sol_small")

# assign these groups into a different column
d$group2 <- d$taxon
d$group2[d$taxon %in% Deep_Hydro] <- "Deep Hydro"
d$group2[d$taxon %in% Liriope] <- "Liriope"
d$group2[d$taxon %in% Shallow_Narco] <- "Shallow Narco"
d$group2[d$taxon %in% Solmaris] <- "Solmaris"
d$group2[d$taxon %in% Other_Hydro] <- "Other Hydro"
d$group2[d$taxon %in% Cydippida] <- "Cydippida"
d$group2[d$taxon %in% Lobata_Thalasso] <- "Lobata_Thalasso"
d$group2[d$taxon %in% Prayidae] <- "Prayidae"

# exclude rare taxa
`%ni%` <- Negate(`%in%`) 
exclude <- c("Charistephane", "Dryodora glandiformis", "Pleurobrachia", "Unknown", "Annatiara", "h10_Pegantha", "h13", "h9_Arctapodema", "r1", "r2", "r3", "r4_Aegina", "r5_Eutonia", "Bolinopsis", "h1")
d <- d[d$group2 %ni% exclude,]

# compute total concentration per group
dg <- ddply(d, ~transect + cast + front + dateTimeB + group2, function(x) {
  tot <- sum(x$concentration)
  timeavg <- mean(x$dateTime)
  return(data.frame(concentration=tot, dateTime=timeavg))
}, .parallel=TRUE)

# convert to wide format
dW <- dcast(dg, dateTimeB~group2, value.var="concentration")

# rename columns
rename <- c("appendicularians" = "Apps", "doliolids" = "Doliolids", "Physonect" = "Physonectae")
dW <- rename(dW, rename)

# set the order
# use the results of the correspondence analysis
# --> CA showed that you can group the taxa into six different groups, 1) velamen, doliolids, liriope, appendicularians, shallow narco, other hydro and solmaris, 2) lobata_thalasso, 3) beroida, 4) cydippida, diphyidae, physonect and sphaeronectes, 5) deep hydro, 6) prayidae
levels <- c("Velamen", "Doliolids", "Liriope", "Apps", "Shallow Narco", "Other Hydro", "Solmaris", "Lobata_Thalasso", "Beroida", "Cydippida", "Diphyidae", "Physonectae", "Sphaeronectes", "Deep Hydro", "Prayidae")

# reorder the columns
dW <- dW[,levels]

# calculate the spearman's correlation coefficient
dcorr <- cor(dW, use="complete.obs", method="spearman")
dc <- melt(dcorr)

# remove the identity values
dc[which(dc$value==1),"value"] <- NA

# save only the lower triangle's rho values (for text later)
dcorr[upper.tri(dcorr,diag=T)] <- NA
dctext <- melt(dcorr)

# plot the heatmap
cortheme <- theme(axis.text.x=element_text(angle=90, vjust=0.5, size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title=element_text(size=26))

ggplot(mapping=aes(x=factor(Var1, levels=levels), y=factor(Var2, levels=levels))) + geom_tile(aes(fill=value), data=dc) + geom_text(aes(label=round(value,2)), data=dctext) + labs(x="", y="", title="", fill="") + scale_fill_gradient2(limits=c(-1,1), low="red", high="blue", na.value="grey90") + theme_bw() + cortheme

# adding in significant values
# try using corrplot package
library(corrplot)
cor.mtest <- function(mat, conf.level = 0.95) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], conf.level = conf.level)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      lowCI.mat[i, j] <- lowCI.mat[j, i] <- tmp$conf.int[1]
      uppCI.mat[i, j] <- uppCI.mat[j, i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(cor(dW, use="complete.obs", method="spearman"), 0.9)
corr <- cor(dW, use="complete.obs", method="spearman")
corrplot(corr, method="square")

corrplot.mixed(corr, upper="color", lower="number")

corrplot(cor(dW, use="complete.obs", method="spearman"), p.mat = res1[[1]], sig.level = 0.1)
corrplot(cor(dW, use="complete.obs", method="spearman"), p.mat = res1[[1]], insig="blank")

# }

## { Construct a correlation heatmap for all taxa without pre-defined groups ---------------------------

# convert to wide format
dW <- dcast(d, dateTimeB~taxon, sum, value.var="concentration")

# rename columns
rename <- c("agel" = "A. elegans", "appendicularians" = "Appendicularians", "doliolids" = "Doliolids", "h11_Haliscera" = "Haliscera sp.2", "h2_Haliscera" = "H. conica", "h3_Cunina" = "Solmaris sp.2", "h5_Liriope" = "L. tetraphylla", "h6_Solmundella" = "S. bitentaculata", "h7_Pegantha" = "Pegantha", "h7_Rhopalonema" = "R. velatum", "h9_Aglaura" = "Aglantha", "Haeckelia beehlri" = "H. beehlri","Hormiphora californiensis" = "H. californiensis", "lemu" = "Diphyidae", "Lilyopsis" = "L. rosea", "Mertensid" = "Mertensiid", "muat" = "M. atlantica", "nabi" = "N. bijuga", "Ocyropsis maculata" = "O. maculata","Solmaris" = "S. rhodoloma", "Thalassocalycidae inconstans" = "T. inconstans", "Velamen" = "V. parallelum")
dW <- rename(dW, rename)
 
# use the results from the ungrouped CA
levels <- c("Pegantha", "Appendicularians", "h15", "S. rhodoloma", "vsh", "S. bitentaculata", "L. tetraphylla", "Doliolids", "V. parallelum", "H. californiensis", "A. elegans", "M. atlantica", "Mertensiid", "N. bijuga", "H. beehlri", "Sphaeronectes", "Beroida", "Larval Lobata", "Prayidae", "T. inconstans", "O. maculata", "L. rosea", "Aglantha", "Diphyidae", "Haliscera sp.2",  "Solmaris sp.2", "R. velatum", "H. conica")

# reorder the columns
dW <- dW[,levels]

# calculate the spearman's correlation coefficient
dcorr <- cor(dW, use="complete.obs", method="spearman")
dc <- melt(dcorr)

# remove the identity values
dc[which(dc$value==1),"value"] <- NA

# save only the lower triangle's rho values (for text later)
dcorr[upper.tri(dcorr,diag=T)] <- NA
dctext <- melt(dcorr)

# plot the heatmap
cortheme <- theme(axis.text.x=element_text(angle=90, vjust=0.5, size=14), axis.text.y=element_text(size=14), legend.text=element_text(size=12), legend.title=element_text(size=12), plot.title=element_text(size=26))

# axis text
axis.text <- expression(italic("Pegantha"), "Appendicularians", "h15", italic("S. rhodoloma"), "vsh", italic("S. bitentaculata"), italic("L. tetraphylla"), "Doliolids", italic("V. parallelum"), italic("H. californiensis"), italic("A. elegans"), italic("M. atlantica"), "Mertensiid", italic("N. bijuga"), italic("H. beehlri"), italic("Sphaeronectes"), "Beroida", "Larval Lobata", "Prayidae", italic("T. inconstans"), italic("O. maculata"), italic("L. rosea"), italic("Aglantha"), "Diphyidae", paste(italic("Haliscera")," sp.2"), paste(italic("Solmaris"), " sp.2"), italic("R. velatum"), italic("H. conica"))

pdf("plots/corr_heatmap/all_withvalues_italics.pdf", width=14, height=13)
p <- ggplot(mapping=aes(x=factor(Var1, levels=levels), y=factor(Var2, levels=levels))) + geom_tile(aes(fill=value), data=dc) + geom_text(aes(label=round(value,2)), size=3, data=dctext) + labs(x="", y="", title="", fill="") + scale_fill_gradient2(limits=c(-1,1), low="red", high="blue", na.value="grey90") + scale_x_discrete(labels=axis.text) + scale_y_discrete(labels=axis.text) + theme_bw() + cortheme
print(p)
dev.off()


# average rho by assemblage
assemA <- c("Pegantha", "h15", "S. rhodoloma", "vsh", "Appendicularians", "S. bitentaculata", "L. tetraphylla", "Doliolids", "V. parallelum", "H. californiensis", "A. elegans")
assemB <- c( "M. atlantica", "Mertensiid", "N. bijuga", "H. beehlri", "Sphaeronectes", "Beroida", "Larval Lobata", "Prayidae")
assemC <- c("T. inconstans", "O. maculata")
assemD <- c("Diphyidae", "L. rosea", "Aglantha", "Solmaris sp.2", "Haliscera sp.2")

avgrho <- data.frame()

dx <- dW[,assemA]
dxc <- cor(dx, use="complete.obs", method="spearman")
diag(dxc) <- NA
avgrho <- rbind(avgrho, data.frame(assemblage="A", avgrho=mean(dxc, na.rm=T)))

dx <- dW[,assemB]
dxc <- cor(dx, use="complete.obs", method="spearman")
diag(xc) <- NA
avgrho <- rbind(avgrho, data.frame(assemblage="B", avgrho=mean(dxc, na.rm=T)))

dx <- dW[,assemC]
dxc <- cor(dx, use="complete.obs", method="spearman")
diag(dxc) <- NA
avgrho <- rbind(avgrho, data.frame(assemblage="C", avgrho=mean(dxc, na.rm=T)))

dx <- dW[,assemD]
dxc <- cor(dx, use="complete.obs", method="spearman")
diag(dxc) <- NA
avgrho <- rbind(avgrho, data.frame(assemblage="D", avgrho=mean(dxc, na.rm=T)))


avgrho

# }
