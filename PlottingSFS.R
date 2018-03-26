setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook")

# Read in frq files
# skip line one and use false for header so it doesn't think that row 2 is a header
IT <- read.table("IT.frq",skip=1,header=F)
WA <- read.table("WA.frq",skip=1,header=F)
NC <- read.table("NC.frq",skip=1,header=F)
head(IT)

# gene name with SNP
# bp site
# alleles
# max 48 (lower is missing data)
# freq

# V1  V2 V3 V4       V5        V6
# 1 OTAU000033-RA  90  2 40 1.000000 0.0000000
# 2 OTAU000033-RA 468  2 44 0.954545 0.0454545
# 3 OTAU000033-RA 510  2 46 0.956522 0.0434783
# 4 OTAU000033-RA 582  2 46 0.869565 0.1304350
# 5 OTAU000033-RA 639  2 48 0.750000 0.2500000
# 6 OTAU000033-RA 717  2 46 0.956522 0.0434783

dim(IT)

# Make a histogram
hist(IT$V6, col="red")
# We want the y axis as a relative frequency, not total number of counts
df <- hist(IT$V6, plot=F)
df$counts <- df$counts/sum(df$counts)
plot(df$counts,type="l")

# Add the other two pops to the plots
dfWA <- hist(WA$V6, plot=F)
dfWA$counts <- dfWA$counts/sum(dfWA$counts)
dfNC <- hist(NC$V6, plot=F)
dfNC$counts <- dfNC$counts/sum(dfNC$counts)

plot(df$counts,type="l",col="red",lwd=1.5,ylim=c(0,0.5))
lines(dfWA$counts,type="l",col="blue",lwd=1.5)
lines(dfNC$counts,type="l",col="green",lwd=1.5)

