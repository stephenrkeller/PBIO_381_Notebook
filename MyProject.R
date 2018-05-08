# My Project
# April 15, 2018
# EAT

setwd("~/Documents/UVM_2018/PBIO381/Project")

# Preliminaries
library("WGCNA")

# Read in the read counts data
beetData = read.table("allcountsdataRN.txt",header=T)
write.csv(beetData, file = "allcountsdataRN.csv",row.names = F)
beetData2 = read.csv("allcountsdataRN.csv")

dim(beetData2)
head(beetData2)
names(beetData2)

bdatExpr0 = as.data.frame(t(beetData[, -c(1:1)]))
dim(bdatExpr0)
names(bdatExpr0) = beetData$ContigName

rownames(bdatExpr0)=names(beetData)[-c(1:1)]
# subset by population
bdatExpr0NC <- subset(bdatExpr0[25:48,])
bdatExpr0WA <- subset(bdatExpr0[49:72,])
bdatExpr0IT <- subset(bdatExpr0[1:24,])
# removing ones with too few counts
bgsg = goodSamplesGenes(bdatExpr0, verbose=3)
bgsgNC = goodSamplesGenes(bdatExpr0NC, verbose=3)
bgsgWA = goodSamplesGenes(bdatExpr0WA, verbose=3)
bgsgIT = goodSamplesGenes(bdatExpr0IT, verbose=3)

################ Make a cluster tree to look for outlier samples

sampleTreeNC = hclust(dist(bdatExpr0NC),method="average")
## plot
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar=c(0,4,2,0))

plot(sampleTreeNC,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

################ Make a cluster tree to look for outlier samples

sampleTreeWA = hclust(dist(bdatExpr0WA),method="average")
## plot
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar=c(0,4,2,0))

plot(sampleTreeWA,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

################ Make a cluster tree to look for outlier samples

sampleTreeIT = hclust(dist(bdatExpr0IT),method="average")
## plot
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar=c(0,4,2,0))

plot(sampleTreeIT,main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab=1.5, cex.axis=1.5, cex.main=2)

######### Read in trait data
# I don't think I will need this information. The file isn't in the correct format so I will try to move forward without that.

########## Network construction and module detection
# create a vector of powers to text
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
powers
# call the network topology function ### this step takes a while to perform!!!

sft  = pickSoftThreshold(bdatExpr0, powerVector = powers, verbose=5)

# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1    0.820  0.274          0.786  4500.0   4330.00   7660
# 2      2    0.504 -0.289          0.615  2080.0   1900.00   4910
# 3      3    0.789 -0.619          0.797  1220.0    983.00   3590
# 4      4    0.801 -0.743          0.808   816.0    566.00   2800
# 5      5    0.782 -0.829          0.794   587.0    342.00   2270
# 6      6    0.785 -0.877          0.808   443.0    217.00   1880
# 7      7    0.771 -0.922          0.797   347.0    147.00   1590
# 8      8    0.762 -0.955          0.797   278.0    101.00   1360
# 9      9    0.756 -0.975          0.789   228.0     74.00   1180
# 10    10    0.745 -0.988          0.773   190.0     55.90   1020
# 11    12    0.892 -0.913          0.886   137.0     32.20    799
# 12    14    0.981 -0.950          0.990   103.0     20.20    693
# 13    16    0.974 -1.050          0.989    80.4     13.50    649
# 14    18    0.964 -1.120          0.969    64.4      9.08    611
# 15    20    0.953 -1.160          0.945    52.8      6.12    578

sftNC = pickSoftThreshold(bdatExpr0NC, powerVector = powers, verbose=5)

# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1    0.489  2.730          0.734  5000.0    5000.0   8350
# 2      2    0.125 -0.134          0.443  2440.0    2190.0   5560
# 3      3    0.785 -0.590          0.752  1480.0    1170.0   4210
# 4      4    0.841 -0.763          0.814  1010.0     690.0   3400
# 5      5    0.840 -0.865          0.826   745.0     435.0   2860
# 6      6    0.818 -0.931          0.819   575.0     287.0   2460
# 7      7    0.768 -1.000          0.779   460.0     196.0   2150
# 8      8    0.739 -1.050          0.767   377.0     142.0   1900
# 9      9    0.732 -1.090          0.774   315.0     109.0   1690
# 10    10    0.729 -1.110          0.776   268.0      85.4   1520
# 11    12    0.715 -1.150          0.772   200.0      54.9   1250
# 12    14    0.694 -1.180          0.740   156.0      38.3   1040
# 13    16    0.683 -1.150          0.688   125.0      27.2    878
# 14    18    0.997 -0.972          0.997   103.0      19.5    749
# 15    20    0.969 -1.010          0.960    86.1      14.1    687

sftWA = pickSoftThreshold(bdatExpr0WA, powerVector = powers, verbose=5)

# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.4880  1.250         0.4700    5730    5680.0   9270
# 2      2   0.0991 -0.152        -0.0823    3220    2580.0   6790
# 3      3   0.4890 -0.431         0.5870    2190    1370.0   5570
# 4      4   0.6800 -0.542         0.7820    1650    1020.0   4790
# 5      5   0.7690 -0.608         0.8260    1300     710.0   4210
# 6      6   0.8110 -0.657         0.8240    1070     508.0   3760
# 7      7   0.8350 -0.698         0.8210     894     386.0   3390
# 8      8   0.8430 -0.731         0.8110     763     295.0   3090
# 9      9   0.8540 -0.766         0.8130     660     230.0   2840
# 10    10   0.8500 -0.795         0.8080     577     182.0   2620
# 11    12   0.8340 -0.844         0.8020     452     129.0   2260
# 12    14   0.8140 -0.898         0.7980     364      90.5   1970
# 13    16   0.7950 -0.943         0.8000     299      64.7   1740
# 14    18   0.7810 -0.982         0.7980     250      50.3   1540
# 15    20   0.7830 -1.010         0.8090     212      38.1   1380

sftIT = pickSoftThreshold(bdatExpr0IT, powerVector = powers, verbose=5)

# Power SFT.R.sq   slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.6540  1.9900          0.787    5680    5800.0   8790
# 2      2   0.0092  0.0242          0.130    3030    2790.0   6170
# 3      3   0.6440 -0.3680          0.631    1960    1810.0   4830
# 4      4   0.7620 -0.5410          0.755    1400    1220.0   3960
# 5      5   0.8480 -0.6430          0.852    1060     854.0   3340
# 6      6   0.8460 -0.7340          0.860     837     619.0   2910
# 7      7   0.8340 -0.8100          0.853     680     455.0   2580
# 8      8   0.8170 -0.8640          0.845     565     338.0   2300
# 9      9   0.8080 -0.9080          0.841     478     254.0   2080
# 10    10   0.7810 -0.9490          0.823     410     196.0   1890
# 11    12   0.7640 -1.0000          0.819     311     128.0   1580
# 12    14   0.7700 -1.0400          0.827     244      84.5   1350
# 13    16   0.7550 -1.0800          0.814     197      58.9   1170
# 14    18   0.7290 -1.1300          0.794     162      42.4   1030
# 15    20   0.7100 -1.1700          0.779     136      31.7    913

####### bdatExpr0
# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

# The power would be 12 for the bdatExpr0. Because it is above 10 should I maybe redo the powers including the odd numbers between 1 and 20? Actually I dont think it would likely make a difference.

####### bdatExpr0NC
# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sftNC$fitIndices[,1], -sign(sftNC$fitIndices[,3])*sftNC$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sftNC$fitIndices[,1], -sign(sftNC$fitIndices[,3])*sftNC$fitIndices[,2], labels = powers, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sftNC$fitIndices[,1], sftNC$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sftNC$fitIndices[,1], sftNC$fitIndices[,5], labels=powers, cex=cex1, col="red")

# The power would be 18 for the bdatExpr0NC. Because it is above 10 should I maybe redo the powers including the odd numbers between 1 and 20? Actually I dont think it would likely make a difference.

####### bdatExpr0WA
# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sftWA$fitIndices[,1], -sign(sftWA$fitIndices[,3])*sftWA$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sftWA$fitIndices[,1], -sign(sftWA$fitIndices[,3])*sftWA$fitIndices[,2], labels = powers, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sftWA$fitIndices[,1], sftWA$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sftWA$fitIndices[,1], sftWA$fitIndices[,5], labels=powers, cex=cex1, col="red")

# There are not any power values that make it above .9 so I should probably re run this with higher power values.

####### bdatExpr0IT
# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sftIT$fitIndices[,1], -sign(sftIT$fitIndices[,3])*sftIT$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sftIT$fitIndices[,1], -sign(sftIT$fitIndices[,3])*sftIT$fitIndices[,2], labels = powers, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sftIT$fitIndices[,1], sftIT$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sftIT$fitIndices[,1], sftIT$fitIndices[,5], labels=powers, cex=cex1, col="red")

# There are not any power values that make it above .9. I don't think running it with higher power values will make a difference though because they already have two peaks 

# Now I will start identifying modules. I will run it for total pop and each of the pops with and without power (or just without power if no power was identified)

# Error: REAL() can only be applied to a 'numeric', not a 'integer'
# ... huh...
# make the columns numerics

bdatExpr0[] <- lapply(bdatExpr0, as.numeric)

netwithpower = blockwiseModules(data.frame(bdatExpr0), power=12, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithpower$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithpower$dendrograms[[1]],mergedColors[netwithpower$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

#
summary(netwithpower)
# Can't correlate trait data so won't do that part.

# Now do bdatExpr0 with a power of 1
netwithoutpower = blockwiseModules(data.frame(bdatExpr0), power=1, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithoutpower$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithoutpower$dendrograms[[1]],mergedColors[netwithoutpower$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# redo with NC with the power of 18
bdatExpr0NC[] <- lapply(bdatExpr0NC, as.numeric)

netwithpowerNC = blockwiseModules(data.frame(bdatExpr0NC), power=18, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithpowerNC$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithpowerNC$dendrograms[[1]],mergedColors[netwithpowerNC$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# redo with NC with the power of 1

netwithoutpowerNC = blockwiseModules(data.frame(bdatExpr0NC), power=18, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithoutpowerNC$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithoutpowerNC$dendrograms[[1]],mergedColors[netwithoutpowerNC$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


# I need to check higher power values for the WA:
powersplus = c(c(1:10), seq(from = 12, to = 30, by = 2))
sftWA = pickSoftThreshold(bdatExpr0WA, powerVector = powersplus, verbose=5)

# Power SFT.R.sq  slope truncated.R.sq mean.k. median.k. max.k.
# 1      1   0.4880  1.250         0.4700    5730    5680.0   9270
# 2      2   0.0991 -0.152        -0.0823    3220    2580.0   6790
# 3      3   0.4890 -0.431         0.5870    2190    1370.0   5570
# 4      4   0.6800 -0.542         0.7820    1650    1020.0   4790
# 5      5   0.7690 -0.608         0.8260    1300     710.0   4210
# 6      6   0.8110 -0.657         0.8240    1070     508.0   3760
# 7      7   0.8350 -0.698         0.8210     894     386.0   3390
# 8      8   0.8430 -0.731         0.8110     763     295.0   3090
# 9      9   0.8540 -0.766         0.8130     660     230.0   2840
# 10    10   0.8500 -0.795         0.8080     577     182.0   2620
# 11    12   0.8340 -0.844         0.8020     452     129.0   2260
# 12    14   0.8140 -0.898         0.7980     364      90.5   1970
# 13    16   0.7950 -0.943         0.8000     299      64.7   1740
# 14    18   0.7810 -0.982         0.7980     250      50.3   1540
# 15    20   0.7830 -1.010         0.8090     212      38.1   1380
# 16    22   0.7780 -1.030         0.8160     182      29.4   1240
# 17    24   0.7690 -1.060         0.8170     158      23.7   1120
# 18    26   0.7660 -1.080         0.8190     138      19.1   1020
# 19    28   0.7600 -1.090         0.8120     122      15.6    924
# 20    30   0.7570 -1.100         0.8060     108      12.7    844

# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sftWA$fitIndices[,1], -sign(sftWA$fitIndices[,3])*sftWA$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sftWA$fitIndices[,1], -sign(sftWA$fitIndices[,3])*sftWA$fitIndices[,2], labels = powersplus, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sftWA$fitIndices[,1], sftWA$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sftWA$fitIndices[,1], sftWA$fitIndices[,5], labels=powersplus, cex=cex1, col="red")

# Still no power hits the 90 line. 9 is the closest though

# redo with WA with the power of 9
bdatExpr0WA[] <- lapply(bdatExpr0WA, as.numeric)

netwithpowerWA = blockwiseModules(data.frame(bdatExpr0WA), power=9, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithpowerWA$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithpowerWA$dendrograms[[1]],mergedColors[netwithpowerWA$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# redo with WA with the power of 1

netwithoutpowerWA = blockwiseModules(data.frame(bdatExpr0WA), power=1, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithoutpowerWA$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithoutpowerWA$dendrograms[[1]],mergedColors[netwithoutpowerWA$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# Rerun all IT for different powers
sftIT = pickSoftThreshold(bdatExpr0IT, powerVector = powersplus, verbose=5)
# plot the results
sizeGrWindow(9,5)
par(mfrow=c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft thresholding power
plot(sftIT$fitIndices[,1], -sign(sftIT$fitIndices[,3])*sftIT$fitIndices[,2], xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type="n", main=paste("Scale independence"))

text(sftIT$fitIndices[,1], -sign(sftIT$fitIndices[,3])*sftIT$fitIndices[,2], labels = powersplus, cex=1, col="red")

abline(h=0.9, col="red")

# mean connectivity as a function of the soft thresholding power

plot(sftIT$fitIndices[,1], sftIT$fitIndices[,5], xlab="Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean Connectivity"))
text(sftIT$fitIndices[,1], sftIT$fitIndices[,5], labels=powersplus, cex=cex1, col="red")

# a power of 30 meets the .9 value


# redo with IT with the power of 
bdatExpr0IT[] <- lapply(bdatExpr0IT, as.numeric)

netwithpowerIT = blockwiseModules(data.frame(bdatExpr0IT), power=30, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithpowerIT$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithpowerIT$dendrograms[[1]],mergedColors[netwithpowerIT$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# redo with IT with the power of 1

netwithoutpowerIT = blockwiseModules(data.frame(bdatExpr0IT), power=1, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0.0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs=TRUE, saveTOMFileBase = "beetleTOM", verbose=3)

sizeGrWindow(12,9)
mergedColors = labels2colors(netwithoutpowerIT$colors)

## plot the dendrogram with the module colors underneath
plotDendroAndColors(netwithoutpowerIT$dendrograms[[1]],mergedColors[netwithoutpowerIT$blockGenes[[1]]],"Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)


###################################################################
# Conclusions:
# I will use data from the pooled populations and a power of 12. This generated 30 modules

library(igraph)
library(WGCNA)
#install.packages("devtools")
library(devtools)
#install_github("jtlovell/limmaDE2")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
biocLite("qvalue")
biocLite("edgeR")
#install.packages("venneuler")
#install.packages("rJava")
library(rJava)
library(limmaDE2)
library(doParallel)
detectCores()

saveRDS(netwithpower, file="netwithpower.rds")


colors <- palette(rainbow(30))
length(colors)
modules2 <- c(0,seq(1:29))
modules2
str(netwithpower)
network <- wgcna2igraph(net = netwithpower, datExpr = as.matrix(bdatExpr0), modules2plot = modules2, colors2plot = colors)

# save network object so i can load it in in the future
#saveRDS(network, file="networkObject.rds")
network <- readRDS(file="networkObject.rds")

plot(network, vertex.label=NA, vertex.size=2,vertex.color="purple",
     edge.color="green") # taking a long time to load...

# It appears to have worked, now I can start calculating centrality measures
# Degree: the number of ties
deg <- degree(network, mode="in")
str(deg)
hist(deg)
plot(deg)
# Another way to get centrality with more information
# res: vertex values(degree, centrality, etc)
# centralization: 
# theoretical_max: miximum score for a graph of this size
# deg2 <- centr_degree(network, mode="in", normalized=T)
# str(deg2)
# hist(deg2[[1]])

# Closeness: centrality based on distance to others in the graph
str(network)
network
# clos <- closeness(network, mode="all") 
# saveRDS(clos, file="closObject.rds")
clos <- readRDS(file = "closObject.rds")

str(clos)
plot(clos)
hist(clos)
# clos2 <- centr_clo(network, mode="all", normalized=T) 
# str(clos2)
# plot(clos2[[1]])
# hist(clos2[[1]])

# Eigenvectors: centrality proportional to the sum of connection centralities
# eig <- eigen_centrality(network, directed=F)
# saveRDS(eig, file="eigObject.rds")
eig <- readRDS(file = "eigObject.rds")

str(eig)
plot(eig[[1]])
hist(eig[[1]])
# eig2 <- centr_eigen(network, directed=T, normalized=T) 
# str(eig2)
# plot(eig2[[1]])
# hist(eig2[[1]])

# Betweenness: centrality based on broker position connecting others
# broker position is a position that is highly connected to groups that are not necessarily connected to each other otherwise
# bet <- betweenness(network, directed=F)
# saveRDS(bet, file="betObject.rds")
bet <- readRDS(file = "betObject.rds")

str(bet)
plot(bet)
hist(bet)
# bet2 <- edge_betweenness(network, directed=F, weights=NA)
# str(bet2)
# plot(bet2)
# hist(bet2)
# bet3 <- centr_betw(network, directed=F, normalized=T)
# str(bet3)
# plot(bet3[[1]])
# hist(bet3[[1]])

# hubs have many outgoing links
# # authorities have many incoming links
# hs <- hub_score(immuno, weights=NA)$vector
# as <- authority_score(immuno, weights=NA)$vector
# par(mfrow=c(1,2))
# plot(immuno, vertex.size=hs*50, main="Hubs")
# plot(immuno, vertex.size=as*30, main="Authorities")

###########################################################################
# Now I have deg, clos, eig, and bet.
setwd("~/Documents/UVM_2018/PBIO381/Project")

library(igraph)
library(WGCNA)
library(devtools)
library(rJava)
library(limmaDE2)

# Read in objects
network <- readRDS(file="networkObject.rds")
deg <- degree(network, mode="in")
clos <- readRDS(file = "closObject.rds")
eig <- readRDS(file = "eigObject.rds")
bet <- readRDS(file = "betObject.rds")

# plot centrality values
par(mfrow=c(2,2))
plot(deg,main = "Degree")
plot(clos, main= "Closeness")
plot(eig[[1]], main = "Eigenvector")
plot(bet, main = "Betweenness")

# create histograms of centrality values
hist(deg,main = "Degree")
hist(clos, main= "Closeness")
hist(eig[[1]], main = "Eigenvector")
hist(bet, main = "Betweenness")
# Eigenvector and betweenness look good because they show distributions with a tail.
# However, betweenness is not a great indicator because it looks at broker positions that don't necessarily tell us too much

# read in putatively deleterious sites based off of snpeff annotations
sites <- read.table(file = "Bad_SNPs_sites.kept.sites", header = T)
head(sites)
class(sites)
dim(sites)
# convert to matrix and put into a vector
sites <- as.matrix(sites)
head(sites)
sitesvec <- vector()
sitesvec <- c(rep(NA, 25013))
head(sitesvec)
length(sitesvec)
counter <- c(rep(1:25013))
head(counter)
for (i in 1:25013) {
  sitesvec[i] <- sites[i,1]
} 
head(sitesvec)
# Extract the tail values from the eigenvector and betweenness
sigeig <- eig[[1]][eig[[1]] >= 0.05]
head(sigeig)
plot(sigeig)
hist(sigeig)
length(sigeig)

sigbet <- bet[bet >= 50000]
head(sigbet)
plot(sigbet)
hist(sigbet)
length(sigbet)

sigeignam <- names(sigeig)
sigbetnam <- names(sigbet)

length(sigeignam)
length(sigbetnam)

sigeignamsnps <- intersect(sites,sigeignam)
length(sigeignamsnps)
sigbetnamsnps <- intersect(sites,sigbetnam)
length(sigbetnamsnps)

sigsnpsitescommon <- intersect(sigeignamsnps,sigbetnamsnps)
length(sigsnpsitescommon)

betgenes <- unique(sigbetnam)
length(betgenes)

# make venn diagram
sitesbet <- cbind(sites,sigeignam)
a <- vennCounts(sitesbet)
install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(x = list(sites, sigeignam),category.names = c("Total Sites","Eigenvector \nSites"), filename = "Sites_Sigeignam.png", main = "Comparison between all putatively deleterious sites and then central eigenvector sites", force.unique = T, fill = c("green", "purple"), cat.pos = c(200, 180))

venn.diagram(x = list(sites, sigbetnam),category.names = c("Total Sites","Betweenness \nSites"), filename = "Sites_Sigbetnam.png", main = "Comparison between all putatively deleterious sites and then central betweenness sites", force.unique = T, fill = c("green", "purple"), cat.pos = c(320, 20), ext.text = F)

venn.diagram(x = list(sigeignamsnps, sigbetnamsnps),category.names = c("Eigenvector \nSites","Betweenness \nSites"), filename = "sigeignamsnps_sigbetnamsnps.png", main = "Comparison between the central eigenvector and betweenness sites", force.unique = T, fill = c("green", "purple"), cat.pos = c(320, 20), ext.text = F)

# make venn diagrams for Jamie
Jamie <- read.table("AlleleAssignmentsBadSNPs.xlsx")


IT200 <- as.matrix(read.table("IT200.txt", header = T))
NC200 <- as.matrix(read.table("NC200.txt", header = T))
WA200 <- as.matrix(read.table("WA200.txt", header = T))
ITQTL <- as.matrix(read.table("ITQTL.txt", header = T))
NCQTL <- as.matrix(read.table("NCQTL.txt", header = T))
WAQTL <- as.matrix(read.table("WAQTL.txt", header = T))

length(IT200)
IT200vec <- vector()
IT200vec <- c(rep(NA, 609))
for (i in 1:609) {
  IT200vec[i] <- IT200[i,1]
} 

length(NC200)
NC200vec <- vector()
NC200vec <- c(rep(NA, 554))
for (i in 1:554) {
  NC200vec[i] <- NC200[i,1]
} 

length(WA200)
WA200vec <- vector()
WA200vec <- c(rep(NA, 603))
for (i in 1:603) {
  WA200vec[i] <- WA200[i,1]
} 

length(ITQTL)
ITQTLvec <- vector()
ITQTLvec <- c(rep(NA, 20))
for (i in 1:20) {
  ITQTLvec[i] <- ITQTL[i,1]
} 

length(NCQTL)
NCQTLvec <- vector()
NCQTLvec <- c(rep(NA, 20))
for (i in 1:20) {
  NCQTLvec[i] <- NCQTL[i,1]
} 

length(WAQTL)
WAQTLvec <- vector()
WAQTLvec <- c(rep(NA, 20))
for (i in 1:20) {
  WAQTLvec[i] <- WAQTL[i,1]
} 
length(ITQTL)
length(NCQTL)
length(WAQTL)
venn.diagram(x = list(IT200, NC200, WA200),category.names = c("IT200","NC200","WA200"), filename = "pop200.png", main = "Comparison across populations for highly expressed genes", force.unique = F, fill = c("green", "purple","orange"), ext.text = F)

venn.diagram(x = list(ITQTL, NCQTL, WAQTL),category.names = c("ITQTL","NCQTL","WAQTL"), filename = "popQTL.png", main = "Comparison across populations for gene expression changes", force.unique = F, fill = c("green", "purple","orange"), ext.text = F)

all200 <- c(IT200vec, NC200vec, WA200vec)
allQTL <- c(ITQTLvec, NCQTLvec, WAQTLvec)

venn.diagram(x = list(all200, allQTL),category.names = c("all200","allQTL"), filename = "all200QTL.png", main = "Comparison of genes identified in each method", force.unique = T, fill = c("green", "purple"), ext.text = F)

# problem, I do not have snp locations just the genes they are in
head(sigbetnamsnps)
df <- as.data.frame(nrow(89), ncol(1))
df <- sigbetnamsnps
df

# get the genes from each population
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

head(WAgenesvec)

ITgenesbet <- intersect(ITgenesvec, sigbetnamsnps)
length(ITgenesbet)

NCgenesbet <- intersect(NCgenesvec, sigbetnamsnps)
length(NCgenesbet)

WAgenesbet <- intersect(WAgenesvec, sigbetnamsnps)
length(WAgenesbet)
length(ITgenesvec)
length(NCgenesvec)
length(WAgenesvec)
ITgeneseig <- intersect(ITgenesvec, sigeignamsnps)
length(ITgeneseig)

NCgeneseig <- intersect(NCgenesvec, sigeignamsnps)
length(NCgeneseig)

WAgeneseig <- intersect(WAgenesvec, sigeignamsnps)
length(WAgeneseig)


# can i get positions out of network? If not I will need to go back to the beginning, use regular expressions to concatinate the position to the contig name and rerun everything

network
vertex_attr(network)

common <- intersect(WA200vec, ITgeneseig)
length(common)



summary(network)
str(network)

ME0vec <- readRDS(file = "ME0vec.rds")
ME0vec <- as.numeric(ME0vec)
bet0 <- betweenness(network, directed=F, vids=ME0vec)

ME0 <- subgraph(network, v=ME0vec)
summary(network)


netwithpower
chooseTopHubInEachModule()

# look at intramodular eigenvalues from the WGCNA or the intermodular eigenvalues from the igraph