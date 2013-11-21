#
#      Create a table of species-physical variable ranges
#
#--------------------------------------------------------------------------

library("plyr")
library("reshape2")
library("ggplot2")
library("gridExtra")

`%ni%` <- Negate(`%in%`) 

# setup R to keep decimal seconds in the times
options("digits.secs"=3)

##{ Read data -------------------------------------------------------------
d <- read.csv("data/all_binned_by_depth.csv", stringsAsFactors=FALSE)
d$dateTime <- as.POSIXct(d$dateTime, tz="America/Los_Angeles")
# }

##{ Process data frame, expand bins ---------------------------------------
# define columns to keep
colToKeep <- c("transect", "cast", "dateTime", "depth", "lat", "long", "temp", "salinity", "swRho", "fluoro", "oxygen", "group", "taxon", "abund", "concentration")

# exclude rare taxa
exclude <- c("Pleurobrachia", "Bolinopsis", "Charistephane", "Dryodora glandiformis", "h1", "r2", "r1", "h13", "h10_Pegantha", "r4_Aegina", "Annatiara", "r3", "r5_Eutonia", "Unknown", "Juvenile Lobata", "h9_Arctapodema", "coor", "foed", "Physonect", "Diphyidae")
d <- d[d$taxon %ni% exclude,]

# set factor order
order <- c("h7_Pegantha", "appendicularians", "h15", "Solmaris", "vsh", "h6_Solmundella", "h5_Liriope", "doliolids", "Velamen", "Hormiphora californiensis", "agel", "muat", "Mertensid", "nabi", "Haeckelia beehlri", "Sphaeronectes",  "Beroida", "Larval Lobata", "Prayidae", "Thalassocalycidae inconstans", "Ocyropsis maculata", "Lilyopsis", "h9_Aglaura",  "lemu", "h11_Haliscera", "h3_Cunina", "h7_Rhopalonema", "h2_Haliscera")

# rename as labels
labels <- c("Pegantha", "Appendicularians", "h15", "S. rhodoloma", "vsh", "S. bitentaculata", "L. tetraphylla", "Doliolids", "V. parallelum", "H. californiensis", "A. elegans", "M. atlantica", "Mertensiid", "N. bijuga", "H. beehlri", "Sphaeronectes", "Beroida", "Larval Lobata", "Prayidae", "T. inconstans", "O. maculata", "L. rosea", "Aglantha", "Diphyidae", "Haliscera sp.2",  "Solmaris sp.2", "R. velatum", "H. conica")


d$taxon <- factor(d$taxon, levels=order, labels=labels)

# separate the appendicularians and the rest of the taxa
dapp <- d[d$taxon=="Appendicularians", colToKeep]
dn <- d[d$taxon!= "Appendicularians", colToKeep]

# copy down rows based on the concentration of organisms per bin
dn <- dn[rep(seq_len(nrow(dn)), round(dn[,"concentration"])), ]
row.names(dn) <- seq(nrow(dn))

# copy down rows for apps
dapp <- dapp[rep(seq_len(nrow(dapp)), round(dapp[,"concentration"]/100)), ]
row.names(dapp) <- seq(nrow(dapp))

dn <- rbind(dn, dapp)

# reverse order for the boxplots because of coord_flip
dn$taxon <- factor(dn$taxon, levels=rev(labels))

# }

##{ Plot boxplots -------------------------------------------------------
# set a plot theme
plottheme <- theme_bw() + theme(axis.text.y=element_blank(), plot.margin = unit(c(1,0,0.5,0), "lines"), legend.position="none", axis.title.y=element_blank(), axis.ticks.y=element_blank()) 

pd <- ggplot(data=dn) + geom_boxplot(aes(x=factor(taxon), y=depth, fill=group), outlier.colour="grey70") + coord_flip() + labs(x="") + scale_fill_brewer("Taxon", palette="Set1") + theme_bw() + theme(legend.position="none", plot.margin=unit(c(1,0,0.5,0), "lines"))

ps <- ggplot(data=dn) + geom_boxplot(aes(x=factor(taxon), y=salinity, fill=group), outlier.colour="grey70") + coord_flip() + scale_fill_brewer("Taxon", palette="Set1") + plottheme

pt <- ggplot(data=dn) + geom_boxplot(aes(x=factor(taxon), y=temp, fill=group), outlier.colour="grey70") + coord_flip() + scale_fill_brewer("Taxon", palette="Set1") + plottheme

po <- ggplot(data=dn) + geom_boxplot(aes(x=factor(taxon), y=oxygen, fill=group), outlier.colour="grey70") + coord_flip() + labs(x="") + scale_fill_brewer("Taxon", palette="Set1") + plottheme 

pf <- ggplot(data=dn) + geom_boxplot(aes(x=factor(taxon), y=fluoro, fill=group), outlier.colour="grey70") + coord_flip() + labs(x="") + scale_fill_brewer("Taxon", palette="Set1") + plottheme + theme(plot.margin=unit(c(1,1,0.5,0), "lines"))

pdf("plots/boxplots_all.pdf", height=8.5, width=11)
grid.arrange(pd, ps, pt, po, pf, nrow=1, widths=c(1.75,1,1,1,1.09))
dev.off()

# }

##{ Frequency distribution table ---------------------------------------
pdf("concentration_range.pdf", width=7, height=10)
print_conc <- adply(order, 1, function(taxon)
{
  df <- data.frame(taxon, 
                   PA=round((length(which(d[d$taxon==taxon,"abund"] != 0))/nrow(d[d$taxon==taxon,]))*100,2),
                   max_Abund=max(d[d$taxon==taxon,"abund"]), 
                   max_Conc=round(max(d[d$taxon==taxon,"concentration"]), 2))
  return(df)
})
# NB: there seems to be a problem with the adply when the d$taxon is an ordered factor. If you get an error when running this code, replace "unique(d$taxon)" in the adply call with "order"

print_conc <- print_conc[,-1]
grid.table(print_conc)
dev.off()

# }