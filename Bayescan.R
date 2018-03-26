# This script analyses my Bayescan output
# March 26, 2018

setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook")

# Read in the snp id information

sites <- read.table("out.kept.sites", header = T)
head(sites)

# Bring in log likelhood data
LL <- read.table("beetles_bayescan10.sel", header = T)
head(LL)

# first column is thinned step in the MCMC model. Each of the Fst values is from each of the populations

# logL       Fst1      Fst2       Fst3
# 50010 -1005119 0.03807107 0.1164735 0.05941002
# 50020 -1005687 0.03784295 0.1151930 0.05975366
# 50030 -1005843 0.03814107 0.1151930 0.06002778
# 50040 -1005310 0.03807489 0.1151930 0.06000667
# 50050 -1005593 0.03807489 0.1144262 0.06016152
# 50060 -1005562 0.03790330 0.1146868 0.06030898

# check to see if we cut off enough of the burnin

plot(LL$logL, type="o", col="gray", ylab="Log-Likelihood")
# if we did not cut out enough, we would see the tail of the burn in increasing up to our "fuzzy caterpillar"

# Plot Fsts
par(mfrow=c(3,1))
plot(LL$Fst1, type="o", col="gray", ylab="Fst1")
plot(LL$Fst2, type="o", col="blue", ylab="Fst2")
plot(LL$Fst3, type="o", col="red", ylab="Fst3")

# read in fst table
Fst <- read.table("beetles_bayescan10_fst.txt")
head(Fst)

# local adaptation has a smal qval, positive alpha and fst that matches...
summary(Fst)
