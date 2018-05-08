# This script analyses my Bayescan output
# March 26, 2018

setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook")

# Read in the snp id information

sites <- read.table("out.kept.sites", header = T)
head(sites)

# Bring in log likelhood data
LL10 <- read.table("bayescan1000010.sel", header = T)
LL20 <- read.table("bayescan1000020.sel", header = T)
LL30 <- read.table("bayescan1000030.sel", header = T)
head(LL10)

# first column is thinned step in the MCMC model. Each of the Fst values is from each of the populations

# logL       Fst1      Fst2       Fst3
# 50010 -1005119 0.03807107 0.1164735 0.05941002
# 50020 -1005687 0.03784295 0.1151930 0.05975366
# 50030 -1005843 0.03814107 0.1151930 0.06002778
# 50040 -1005310 0.03807489 0.1151930 0.06000667
# 50050 -1005593 0.03807489 0.1144262 0.06016152
# 50060 -1005562 0.03790330 0.1146868 0.06030898

# check to see if we cut off enough of the burnin
par(mfrow=c(3,1))
plot(LL10$logL, type="o", col="gray", ylab="Log-Likelihood",xlim = c(0,500))
plot(LL20$logL, type="o", col="blue", ylab="Log-Likelihood",xlim = c(0,500))
plot(LL30$logL, type="o", col="red", ylab="Log-Likelihood",xlim = c(0,500))


# if we did not cut out enough, we would see the tail of the burn in increasing up to our "fuzzy caterpillar"

# Plot Fsts
par(mfrow=c(3,1))
plot(LL10$Fst1, type="o", col="gray", ylab="Fst1")
plot(LL10$Fst2, type="o", col="blue", ylab="Fst2")
plot(LL10$Fst3, type="o", col="red", ylab="Fst3")

par(mfrow=c(3,1))
plot(LL20$Fst1, type="o", col="gray", ylab="Fst1")
plot(LL20$Fst2, type="o", col="blue", ylab="Fst2")
plot(LL20$Fst3, type="o", col="red", ylab="Fst3")

par(mfrow=c(3,1))
plot(LL30$Fst1, type="o", col="gray", ylab="Fst1")
plot(LL30$Fst2, type="o", col="blue", ylab="Fst2")
plot(LL30$Fst3, type="o", col="red", ylab="Fst3")

# From now on we are moving forward with the thinning value of 10
# read in fst table
Fst <- read.table("bayescan1000010_fst.txt")
head(Fst)

# local adaptation has a smal qval, positive alpha and fst that matches...
# log 10 Posterior odds --> higher indicate selection, lower mean nuetrality (sometimes called a bayes value)
# alpha is our measure or selection (near 0 nuetral, greater than 0 mean local adaptation selection, negative means stabilizing selction)
# Fst follows beta parameter when alpha is low
# significant values have low q values (similar to p values)
summary(Fst)

# Bring in out.kept.sites file that has positions of our things in the same order
sites <- read.table("out.kept.sites", header = T)

Fst_sites <- cbind(sites,Fst)
head(Fst_sites)
dim(Fst_sites)
dim(sites)
dim(Fst)

# take the q values and log 10 them to spread out the values a little bit for better graphing
Fst_sites$log10q <- -log10(Fst_sites$qval)
head(Fst_sites)

# Now we want to set our false discovery rate, proportion of significant results that are expected to be falst positives
FDR <- 0.01
candsnps_pr10000 <- Fst_sites[which(Fst_sites$log10q>(-log10(FDR))),]
dim(candsnps_pr10000)
# [1] 2 8
# only two significant ;(

# View significant locations
table(droplevels(candsnps_pr10000$CHROM))
dev.off()
plot(Fst_sites$fst, Fst_sites$alpha, xlab="Fst",ylab="alpha")
points(candsnps_pr10000$fst,candsnps_pr10000$alpha,col="red",pch=16)

# write out the candidate snps to annotate them on the server
write.table(candsnps_pr10000$CHROM,"bayescan_pr10000_candsnps.txt",quote = F,row.names = F,col.names = F)
