setwd("~/Documents/UVM_2018/PBIO381/Project")

# Preliminaries
library(WGCNA)
library(igraph)

# read in network object, eigenvalues for the network and the bad snps list
network <- readRDS(file = "networkObject.rds")
eig <- readRDS(file = "eigObject.rds")
badsnpsites <- read.table(file = "Bad_SNPs_sites.kept.sites", header = T)
# put chrom into vector
badsnpvec <- as.vector(badsnpsites$CHROM)
# use intersect to pull eig values for contigs with deleterious snps
eigdelnam <- intersect(badsnpvec, names(eig[[1]]))
# get eig values for those
eigdel <- eig[[1]][eigdelnam]
length(eigdel)
# plot histogram
hist(eigdel)
plot(eigdel)
# make a vector with the eig values for all the genes that are not deleterious
eignotdelnam <- setdiff(names(eig[[1]]),badsnpvec)
eignotdel <- eig[[1]][eignotdelnam]
length(eignotdel)
hist(eignotdel)
plot(eignotdel)

# now lets try doing this per population
# pull in genes with snps from each population
ITgenes <- as.matrix(read.table("ITgenes", header=T))
NCgenes <- as.matrix(read.table("NCgenes", header=T))
WAgenes <- as.matrix(read.table("WAgenes", header=T))

length(ITgenes)
ITgenesvec <- vector()
ITgenesvec <- c(rep(NA, 117813))
for (i in 1:117813) {
  ITgenesvec[i] <- ITgenes[i,1]
} 

length(NCgenes)
NCgenesvec <- vector()
NCgenesvec <- c(rep(NA, 98579))
for (i in 1:98579) {
  NCgenesvec[i] <- NCgenes[i,1]
} 

length(WAgenes)
WAgenesvec <- vector()
WAgenesvec <- c(rep(NA, 108029))
for (i in 1:108029) {
  WAgenesvec[i] <- WAgenes[i,1]
} 

# remove duplicates that arise from multiple snps per contig
ITunigenes <- unique(ITgenesvec)
NCunigenes <- unique(NCgenesvec)
WAunigenes <- unique(WAgenesvec)
length(ITunigenes)

# Now pull out genes with deleterious snps
ITdelnam <- intersect(badsnpvec, ITunigenes)
length(ITdelnam)
NCdelnam <- intersect(badsnpvec, NCunigenes)
WAdelnam <- intersect(badsnpvec, WAunigenes)

# Pull out eig values for genes with deleterious snps per population
ITeigdelnam <- intersect(ITdelnam, names(eig[[1]]))
length(ITeigdelnam)
NCeigdelnam <- intersect(NCdelnam, names(eig[[1]]))
length(NCeigdelnam)
WAeigdelnam <- intersect(WAdelnam, names(eig[[1]]))
length(WAeigdelnam)
# get the eig values for those deleterious contigs per pop
ITeigdel <- eig[[1]][ITeigdelnam]
NCeigdel <- eig[[1]][NCeigdelnam]
WAeigdel <- eig[[1]][WAeigdelnam]
length(ITeigdel)
# Hist and plot these!!!
plot(ITeigdel)
plot(NCeigdel)
plot(WAeigdel)
hist(ITeigdel)
hist(NCeigdel)
hist(WAeigdel)
# now do the same for the non deleterious contigs
ITeignotdelnam <- setdiff(names(eig[[1]]),ITeigdelnam)
ITeignotdel <- eig[[1]][ITeignotdelnam]
NCeignotdelnam <- setdiff(names(eig[[1]]),NCeigdelnam)
NCeignotdel <- eig[[1]][NCeignotdelnam]
WAeignotdelnam <- setdiff(names(eig[[1]]),WAeigdelnam)
WAeignotdel <- eig[[1]][WAeignotdelnam]
# plto and hist!
plot(ITeignotdel)
plot(NCeignotdel)
plot(WAeignotdel)
hist(ITeignotdel)
hist(NCeignotdel)
hist(WAeignotdel)


ks.test(ITeigdel,ITeignotdel)
# lets make some venn diagrams to see how these deleterious genes overlap across the populations

venn.diagram(x = list(names(ITeigdel), names(NCeigdel), names(WAeigdel)),category.names = c("ITeigdel","NCeigdel","WAeigdel"), filename = "popeigdel.png", main = "Comparison across populations for eigdels", force.unique = T, fill = c("green", "purple","orange"), ext.text = F)

# Now let's look at the genes that have high levels of eigenvalues across the three populations
sigITeigdel <- ITeigdel[ITeigdel >= 0.95]
sigNCeigdel <- NCeigdel[NCeigdel >= 0.95]
sigWAeigdel <- WAeigdel[WAeigdel >= 0.95]
# plot and hist
plot(sigITeigdel)
plot(sigNCeigdel)
plot(sigWAeigdel)
hist(sigITeigdel)
hist(sigNCeigdel)
hist(sigWAeigdel)
length(sigITeigdel)
length(sigNCeigdel)
length(sigWAeigdel)

# venn diagram
venn.diagram(x = list(names(sigITeigdel), names(sigNCeigdel), names(sigWAeigdel)),category.names = c("sigITeigdel","sigNCeigdel","sigWAeigdel"), filename = "popsigeigdel.png", main = "Comparison across populations for sigeigdels", force.unique = T, fill = c("green", "purple","orange"), ext.text = F)

# lets get all of the gene names for these. From the venn diagram I can see that IT has all of the genes but one so I will get that
one <- setdiff(names(sigNCeigdel),names(sigITeigdel))
length(one)
one <- sigNCeigdel[one]
one
allsigeig <- c(sigITeigdel, one)
length(sigITeigdel)
length(allsigeig)
head(allsigeig)
allsignam <- names(allsigeig)
allsignamdf <- matrix(nrow = 196, ncol = 1)
allsignamdf[,1] <- allsignam
allsignamdf <- as.data.frame(allsignamdf)
head(allsignamdf)
write.csv(allsignamdf, file = "allsignam")
