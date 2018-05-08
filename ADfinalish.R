setwd("~/Documents/UVM_2018/PBIO381/Project")

library(dplyr)
#Ancestral Derived Here I come
ADpop <- read.csv(file = "AncestralDerivedAlleles.csv", header = T)
Consensus <- read.csv(file = "AncestralDerivedAllelesNew.csv", header = T)
head(Consensus)
dim(Consensus)
class(ADpop)

# Split into pops and consensus
pops <- subset(ADpop, select=c(CHPS,ALTIT,CHPSNC,ALTNC,CHPSWA,ALTWA))
head(pops)
dim(pops)
dim(AD)
# remove extras

pops <- merge(Consensus, pop, by="CHPS")
dim(pops)
head(pops)

# pull out columns i need
df <- subset(pops, select = c(CHPS,Otau.Ref.Allele,Consensus.Ancestral.Allele,ALTIT,ALTNC,ALTWA))
head(df)

df2 <- rbind(subset(df, df$Consensus.Ancestral.Allele == "A"),subset(df, df$Consensus.Ancestral.Allele == "T"),subset(df, df$Consensus.Ancestral.Allele == "C"),subset(df, df$Consensus.Ancestral.Allele == "G"))
dim(df2)
head(df2)
class(df2)
df2 <- df2[order(df2$CHPS),]
head(df2, n=50)

# Now I am ready to begin filtering through to assign ancestral or derived per population for each snp.
dim(df2)
df2 <- cbind(df2, IT=c(rep(NA,126329)))
# accidentally added a column i didn't want
df2 <- df2[,-7]
head(df2)
# add the other two
df2 <- cbind(df2, NC=c(rep(NA,126329)))
df2 <- cbind(df2, WA=c(rep(NA,126329)))
head(df2)

# now start filtering for derived and ancestral
# if the reference is the same as the consensus and pop is different from the consensus --> derived
# if the reference is different from the consensus and the pop is the same as the consensus --> ancestral
# if the reference is different from the consensus and the pop is different from the consensus --> derived
# if the reference is the same as the consensus and the pop is the same as the consensus --> ancestral
# Let's start with IT
for (i in 1:126329) {
  if (df2[i,2] == df2[i,3] | df2[i,4] != df2[i,3]) {
    df2[i,7] <- "D"
  }
  else {
    if (df2[i,2] != df2[i,3] | df2[i,4] == df2[i,3]) {
      df2[i,7] <- "A"
    }
    else {
      if (df2[i,2] != df2[i,3] | df2[i,4] != df2[i,3]) {
        df2[i,7] <- "D"
      }
      else {
        if (df2[i,2] == df2[i,3] | df2[i,4] == df2[i,3])
          df2[i,7] <- "A"
      }
      }
    }
    }


# got an error in different levels of factors so I will try to change snp columns to characters
class(df2$ALTIT)
df2$ALTIT <- as.character(df2$ALTIT)
df2$ALTNC <- as.character(df2$ALTNC)
df2$ALTWA <- as.character(df2$ALTWA)
# try to make it a matrix to see if that works because the above did not
df2 <- as.matrix(df2)
head(df2, n=50)
str(df2)
# OMG it worked!!! I just had forgotten to put 1:126329 instead of just 126329
# Now I will do that for NC and WA
# NC
for (i in 1:126329) {
  if (df2[i,2] == df2[i,3] | df2[i,5] != df2[i,3]) {
    df2[i,8] <- "D"
  }
  else {
    if (df2[i,2] != df2[i,3] | df2[i,5] == df2[i,3]) {
      df2[i,8] <- "A"
    }
    else {
      if (df2[i,2] != df2[i,3] | df2[i,5] != df2[i,3]) {
        df2[i,8] <- "D"
      }
      else {
        if (df2[i,2] == df2[i,3] | df2[i,5] == df2[i,3])
          df2[i,8] <- "A"
      }
    }
  }
}
head(df2)
# Now WA
for (i in 1:126329) {
  if (df2[i,2] == df2[i,3] | df2[i,6] != df2[i,3]) {
    df2[i,9] <- "D"
  }
  else {
    if (df2[i,2] != df2[i,3] | df2[i,6] == df2[i,3]) {
      df2[i,9] <- "A"
    }
    else {
      if (df2[i,2] != df2[i,3] | df2[i,6] != df2[i,3]) {
        df2[i,9] <- "D"
      }
      else {
        if (df2[i,2] == df2[i,3] | df2[i,6] == df2[i,3])
          df2[i,9] <- "A"
      }
    }
  }
}
head(df2, n=150)

# i am concerned that there are no differences in in snps between pops so I will try subsetting to see if the dimensions are different
df2 <- as.data.frame(df2)
head(df2)
saveRDS(df2, file = "DerivedAncestralobject.rds")
ITD <- subset(df2, df2$IT == "D", select = c(CHPS,IT))
head(ITD)
dim(ITD)
NCD <- subset(df2, df2$NC == "D", select = c(CHPS,NC))
head(NCD)
dim(NCD)
WAD <- subset(df2, df2$WA == "D", select = c(CHPS,WA))
head(WAD)
dim(WAD)
# hmmmm that is suspicious but I will move forward and ask about it tomorrow when I meet with Jamie and Steve and Melissa
# alright game plan time. I have a list of my central genes. I have a list of genes and positions with frequences. I have a list of genes and positions that are derived. My final list needs to be a list of central genes that are derived with frequencies above 0.85 and below 0.15 for ones that have gone to fixation and ones that are buing selected against for each population.
# first I need to compare this list of snps to the list of deleterious snps of high and moderate severity.
# read in bad snps
badsnpsites <- read.table(file = "Bad_SNPs_sites.kept.sites", header = T)
head(badsnpsites)
dim(badsnpsites)
badsnpsitesvec <- as.vector(paste(badsnpsites$CHROM,badsnpsites$POS,sep = "_"))
length(badsnpsitesvec)
derITsites <- as.vector(ITD$CHPS)
derbadsnpsitesIT <- intersect(badsnpsitesvec,derITsites)
length(derbadsnpsitesIT)

derNCsites <- as.vector(NCD$CHPS)
derbadsnpsitesNC <- intersect(badsnpsitesvec,derNCsites)
length(derbadsnpsitesNC)

derWAsites <- as.vector(WAD$CHPS)
derbadsnpsitesWA <- intersect(badsnpsitesvec,derWAsites)
length(derbadsnpsitesWA)

# now that i have pop lists for derived and deleterious, I need to split them by frequencies
ITfrq <- read.table(file = "IT.frq", header = T)
head(ITfrq)
# The positions read in as row names and I need them to be a column
ITfrq$CHPS <- rownames(ITfrq)
head(ITfrq)
IThighfrq <- subset(ITfrq, X.FREQ. > 0.15, select = c(CHPS,X.FREQ.))
head(IThighfrq)
# now the low freq (less than 0.15
ITlowfrq <- subset(ITfrq, X.FREQ. <= 0.15, select = c(CHPS, X.FREQ.))
head(ITlowfrq)

# Repeat for NC and WA
NCfrq <- read.table(file = "NC.frq", header = T)
head(NCfrq)
# The positions read in as row names and I need them to be a column
NCfrq$CHPS <- NCfrq$CHROM
head(NCfrq)
NChighfrq <- subset(NCfrq, X.FREQ. > 0.15, select = c(CHPS,X.FREQ.))
head(NChighfrq)
# now the low freq (less than 0.15
NClowfrq <- subset(NCfrq, X.FREQ. <= 0.15, select = c(CHPS, X.FREQ.))
head(NClowfrq)

#WA
WAfrq <- read.table(file = "WA.frq", header = T)
head(WAfrq)
# The positions read in as row names and I need them to be a column
WAfrq$CHPS <- rownames(WAfrq)
head(WAfrq)
WAhighfrq <- subset(WAfrq, X.FREQ. > 0.15, select = c(CHPS,X.FREQ.))
head(WAhighfrq)
# now the low freq (less than 0.15
WAlowfrq <- subset(WAfrq, X.FREQ. <= 0.15, select = c(CHPS, X.FREQ.))
head(WAlowfrq)


dim(IThighfrq)
dim(ITlowfrq)
dim(NChighfrq)
dim(NClowfrq)
dim(WAhighfrq)
dim(WAlowfrq)

# now I need to make vectors of the positions
IThighfrqsites <- as.vector(IThighfrq$CHPS)
NChighfrqsites <- as.vector(NChighfrq$CHPS)
WAhighfrqsites <- as.vector(WAhighfrq$CHPS)

ITlowfrqsites <- as.vector(ITlowfrq$CHPS)
NClowfrqsites <- as.vector(NClowfrq$CHPS)
WAlowfrqsites <- as.vector(WAlowfrq$CHPS)


# now I want to compare my lists of derived and deleterious to the frequency positions
IThighfrqderbad <- intersect(IThighfrqsites, derbadsnpsitesIT)
length(IThighfrqderbad)
NChighfrqderbad <- intersect(NChighfrqsites, derbadsnpsitesNC)
length(NChighfrqderbad)
WAhighfrqderbad <- intersect(WAhighfrqsites, derbadsnpsitesWA)
length(WAhighfrqderbad)

ITlowfrqderbad <- intersect(ITlowfrqsites, derbadsnpsitesIT)
length(ITlowfrqderbad)
NClowfrqderbad <- intersect(NClowfrqsites, derbadsnpsitesNC)
length(NClowfrqderbad)
WAlowfrqderbad <- intersect(WAlowfrqsites, derbadsnpsitesWA)
length(WAlowfrqderbad)

# wow... so that is not what was expected. We expect more high freq deleterious derived snps in NC than WA which had one more than IT and for low freq deleterious snps, IT had the most, then WA, then NC. and we got the opposite

# because these numbers are so small though, I am nervous that none of my central genes will even map on...
# but in order to do that I need to split the contig from the snp because my centrality picked out only contig names not snp positions as well.
# to do this I need to put all of these vectors into a dataframe so I can split the names

IThighfrqderbaddf <- as.data.frame(matrix(data = IThighfrqderbad, ncol = 1))
head(IThighfrqderbaddf)
library(stringr)
IThighfrqderbadcontigdf <- str_split_fixed(IThighfrqderbaddf$V1, "_", 2)
class(IThighfrqderbadcontigdf)
IThighfrqderbadcontig <- as.vector(IThighfrqderbadcontigdf[,1])
IThighfrqderbadcontig

# It worked! now do that for NC and WA and then for all three pops for the low frequencies
NChighfrqderbaddf <- as.data.frame(matrix(data = NChighfrqderbad, ncol = 1))
head(NChighfrqderbaddf)
NChighfrqderbadcontigdf <- str_split_fixed(NChighfrqderbaddf$V1, "_", 2)
head(NChighfrqderbadcontigdf)
NChighfrqderbadcontig <- as.vector(NChighfrqderbadcontigdf[,1])
head(NChighfrqderbadcontig)

# high frq WA
WAhighfrqderbaddf <- as.data.frame(matrix(data = WAhighfrqderbad, ncol = 1))
head(WAhighfrqderbaddf)
WAhighfrqderbadcontigdf <- str_split_fixed(WAhighfrqderbaddf$V1, "_", 2)
head(WAhighfrqderbadcontigdf)
WAhighfrqderbadcontig <- as.vector(WAhighfrqderbadcontigdf[,1])
head(WAhighfrqderbadcontig)

# low frq for all!
# IT
ITlowfrqderbaddf <- as.data.frame(matrix(data = ITlowfrqderbad, ncol = 1))
head(ITlowfrqderbaddf)
ITlowfrqderbadcontigdf <- str_split_fixed(ITlowfrqderbaddf$V1, "_", 2)
head(ITlowfrqderbadcontigdf)
ITlowfrqderbadcontig <- as.vector(ITlowfrqderbadcontigdf[,1])
head(ITlowfrqderbadcontig)

# low freq NC
NClowfrqderbaddf <- as.data.frame(matrix(data = NClowfrqderbad, ncol = 1))
head(NClowfrqderbaddf)
NClowfrqderbadcontigdf <- str_split_fixed(NClowfrqderbaddf$V1, "_", 2)
head(NClowfrqderbadcontigdf)
NClowfrqderbadcontig <- as.vector(NClowfrqderbadcontigdf[,1])
head(NClowfrqderbadcontig)

# low frq WA
WAlowfrqderbaddf <- as.data.frame(matrix(data = WAlowfrqderbad, ncol = 1))
head(WAlowfrqderbaddf)
WAlowfrqderbadcontigdf <- str_split_fixed(WAlowfrqderbaddf$V1, "_", 2)
head(WAlowfrqderbadcontigdf)
WAlowfrqderbadcontig <- as.vector(WAlowfrqderbadcontigdf[,1])
head(WAlowfrqderbadcontig)

# great, so now I need to compare with my centrality genes for each population
eig <- readRDS(file = "eigObject.rds")
str(eig)
length(eig[[1]])
sigeig <- eig[[1]][eig[[1]] >= 0.95]
head(sigeig)
length(sigeig)
sigeignam <- names(sigeig)
head(sigeignam)

# lets see if there is overlap!
IThighfrqderbadeig <- intersect(IThighfrqderbadcontig,sigeignam)
length(IThighfrqderbadeig)
NChighfrqderbadeig <- intersect(NChighfrqderbadcontig,sigeignam)
length(NChighfrqderbadeig)
WAhighfrqderbadeig <- intersect(WAhighfrqderbadcontig,sigeignam)
length(WAhighfrqderbadeig)

ITlowfrqderbadeig <- intersect(ITlowfrqderbadcontig,sigeignam)
length(ITlowfrqderbadeig)
NClowfrqderbadeig <- intersect(NClowfrqderbadcontig,sigeignam)
length(NClowfrqderbadeig)
WAlowfrqderbadeig <- intersect(WAlowfrqderbadcontig,sigeignam)
length(WAlowfrqderbadeig)
# There is no overlap among any of those of high freq but the low freq doesnt match our ideas regarding demographic history... ;( maybe I should decrease the stringency... taking above 0.95 gives me 330 while if I were to take the top 5% it would give me the top 469.35...

# How do I take the top 470 contigs... 
# I could order the vector from smallest to largest and then take a tail of size 470 and put that into a new vector...
orderedeig <- sort(eig[[1]], decreasing = F)
head(orderedeig)
tail(orderedeig)
sigeig2 <- tail(orderedeig,n=470)
tail(sigeig2)
head(sigeig2)
hist(sigeig2)
sigeignam2 <- names(sigeig2)
head(sigeignam2)
length(sigeignam2)

# repeat above but with sigeig2
IThighfrqderbadeig2 <- intersect(IThighfrqderbadcontig,sigeignam2)
length(IThighfrqderbadeig2)
NChighfrqderbadeig2 <- intersect(NChighfrqderbadcontig,sigeignam2)
length(NChighfrqderbadeig2)
WAhighfrqderbadeig2 <- intersect(WAhighfrqderbadcontig,sigeignam2)
length(WAhighfrqderbadeig2)

ITlowfrqderbadeig2 <- intersect(ITlowfrqderbadcontig,sigeignam2)
length(ITlowfrqderbadeig2)
NClowfrqderbadeig2 <- intersect(NClowfrqderbadcontig,sigeignam2)
length(NClowfrqderbadeig2)
WAlowfrqderbadeig2 <- intersect(WAlowfrqderbadcontig,sigeignam2)
length(WAlowfrqderbadeig2)

length(WAlowfrqderbadcontig)

# still didnt get any.
# actually it turns out that I needed to make the high freq ones be above 0.15 not 0.85 and now we get results


# time to do a chi squared test
# make dataframe
# I did this outside of R

# Time to read in Lucy's list of genes under selective sweeps based on their analyses
candsIT <- readRDS(file = "candsITsweed.rds")
head(candsIT)
candsNC <- readRDS(file = "candsNCsweed.rds")
head(candsNC)
candsWA <- readRDS(file = "candsWAsweed.rds")
head(candsWA)
# So the position is coming in as the row names but that is okay because I just need the contig names in order to overlap with my central genes. These are genes generated from SweeD

candsITvec <- as.vector(candsIT$CHROM)
head(candsITvec)
candsNCvec <- as.vector(candsNC$CHROM)
head(candsNCvec)
candsWAvec <- as.vector(candsWA$CHROM)
head(candsWAvec)

# Let's compare them against my high and low freq lists of central, derived, deleterious contigs
ITsweepcentralhf <- intersect(candsITvec, IThighfrqderbadeig2)
ITsweepcentralhf
# 1
NCsweepcentralhf <- intersect(candsNCvec, NChighfrqderbadeig2)
NCsweepcentralhf
# 0
WAsweepcentralhf <- intersect(candsWAvec, WAhighfrqderbadeig2)
WAsweepcentralhf
# 0

ITsweepcentrallf <- intersect(candsITvec, ITlowfrqderbadeig2)
ITsweepcentrallf
# 3
NCsweepcentrallf <- intersect(candsNCvec, NClowfrqderbadeig2)
NCsweepcentrallf
# 0
WAsweepcentrallf <- intersect(candsWAvec, WAlowfrqderbadeig2)
WAsweepcentrallf
# 0

# interesting... so only the selective sweeps identified in IT show overlap with alleles expected to be under BGS based on centrality

# Fst outliers
Fst <- read.table(file = "FSToutliers_Bayscan.txt", header = T)
head(Fst)
Fstvec
# XTX outliers
XTX <- read.table(file = "XtXoutliers_Bayenv2.txt", header = T)
head(XTX)
XTXvec
# They need to become vectors of contigs
Fstvec <- as.vector(Fst$CHROM)
head(Fstvec)
XTXvec <- as.vector(XTX$V1)
head(XTXvec)

# compare lists to get overlap per pop with my central deleterious derived high and low freq lists
ITFsthf <- intersect(IThighfrqderbadeig2, Fstvec)
ITFsthf
ITXTXhf <- intersect(IThighfrqderbadeig2, XTXvec)
ITXTXhf

NCFsthf <- intersect(NChighfrqderbadeig2, Fstvec)
NCFsthf
NCXTXhf <- intersect(NChighfrqderbadeig2, XTXvec)
NCXTXhf

WAFsthf <- intersect(WAhighfrqderbadeig2, Fstvec)
WAFsthf
WAXTXhf <- intersect(WAhighfrqderbadeig2, XTXvec)
WAXTXhf

# low freq
ITFstlf <- intersect(ITlowfrqderbadeig2, Fstvec)
ITFstlf
ITXTXlf <- intersect(ITlowfrqderbadeig2, XTXvec)
ITXTXlf

NCFstlf <- intersect(NClowfrqderbadeig2, Fstvec)
NCFstlf
NCXTXlf <- intersect(NClowfrqderbadeig2, XTXvec)
NCXTXlf

WAFstlf <- intersect(WAlowfrqderbadeig2, Fstvec)
WAFstlf
WAXTXlf <- intersect(WAlowfrqderbadeig2, XTXvec)
WAXTXlf

# seems as though the overlap doesn't usually change per population (is typically the same) but also not all the places identified in selective sweeps

# I also need to go back and get values for chi square analyses
# First intersect I did was between badsnpsites and derived sites
# This means that first I need to intersect the neutral sites with the derived sites... or do I need to intersect the bad snps with the ancestral sites...
# Let's subset the A from my ancestral derived table
ITA <- subset(df2, df2$IT == "A", select = c(CHPS,IT))
head(ITA)
dim(ITA)
NCA <- subset(df2, df2$NC == "A", select = c(CHPS,NC))
head(NCA)
dim(NCA)
WAA <- subset(df2, df2$WA == "A", select = c(CHPS,WA))
head(WAA)
dim(WAA)
# put the positions into vectors
ancITsites <- as.vector(ITA$CHPS)
ancNCsites <- as.vector(NCA$CHPS)
ancWAsites <- as.vector(WAA$CHPS)
# lets see how many deleterious ancestral snps there are
ancbadsnpsitesIT <- intersect(badsnpsitesvec,ancITsites)
ancbadsnpsitesNC <- intersect(badsnpsitesvec,ancNCsites)
ancbadsnpsitesWA <- intersect(badsnpsitesvec,ancWAsites)
length(ancbadsnpsitesIT)
length(ancbadsnpsitesNC)
length(ancbadsnpsitesWA)
length(derbadsnpsitesIT)
length(derbadsnpsitesNC)
length(derbadsnpsitesWA)
# or maybe what I want is how many derived snp sites aren't deleterious... this shows that of all snps, most are deleterious
length(derbadsnpsitesIT)
ITDnotdel <- length(derITsites) - length(derbadsnpsitesIT)
ITDnotdel
length(derbadsnpsitesNC)
NCDnotdel <- length(derNCsites) - length(derbadsnpsitesNC)
NCDnotdel
length(derbadsnpsitesWA)
WADnotdel <- length(derWAsites) - length(derbadsnpsitesWA)
WADnotdel

# this is all garbage... what am I really trying to understand...
# I have "derived snps" which is like my raw snp calls but more realistic. Then I map my deleterious snps on. When I do this I cannot tell which pop has more or less deleterious snps because freq is the only way to weed out ones that are at a freq of 1 for the ref allele. So really, I want to do a chi square for high and low freq with my central contigs mapped on to say, of the derived, deleterious snps at high freq, is there a sig difference between the number of these that are central vs non central across the populations and then again for low freq. 

# Then when I map on the genes identified by selective sweeps... of the select

# Derived neutral HF and LF versus derived deleterious HF and LF
# another version of that with central

# Analyses left to do:

# MWU test on difference in overlap between deleterious and my central and neutral and my central

totalsweepcontigs <- c(Fstvec,XTXvec,candsITvec,candsNCvec,candsWAvec)
unitotalsweepcontigs <- unique(totalsweepcontigs) 
length(unitotalsweepcontigs)
ITtotalsweep <- c(Fstvec,XTXvec,candsITvec)
uniITtotalsweep <- unique(ITtotalsweep)
length(uniITtotalsweep)

NCtotalsweep <- c(Fstvec,XTXvec,candsNCvec)
uniNCtotalsweep <- unique(NCtotalsweep)
length(uniNCtotalsweep)

WAtotalsweep <- c(Fstvec,XTXvec,candsWAvec)
uniWAtotalsweep <- unique(WAtotalsweep)
length(uniWAtotalsweep)

# neutral
neutral <- read.table(file = "Neutral_SNPs_sites.kept.sites.txt", header = T)
head(neutral)
neutralsnpsitesvec <- as.vector(paste(neutral$CHROM,neutral$POS,sep = "_"))
head(neutralsnpsitesvec)
neuderIT <- intersect(derITsites, neutralsnpsitesvec)
neuderNC <- intersect(derNCsites, neutralsnpsitesvec)
neuderWA <- intersect(derWAsites, neutralsnpsitesvec)
neuderhfIT <- intersect(neuderIT, IThighfrqsites)
neuderhfNC <- intersect(neuderNC, NChighfrqsites)
neuderhfWA <- intersect(neuderWA, WAhighfrqsites)
neuderlfIT <- intersect(neuderIT, ITlowfrqsites)
neuderlfNC <- intersect(neuderNC, NClowfrqsites)
neuderlfWA <- intersect(neuderWA, WAlowfrqsites)

length(neuderhfIT)
length(neuderhfNC)
length(neuderhfWA)
length(neuderlfIT)
length(neuderlfNC)
length(neuderlfWA)

neuderhfITdf <- as.data.frame(matrix(data = neuderhfIT, ncol = 1))
neuderhfNCdf <- as.data.frame(matrix(data = neuderhfNC, ncol = 1))
neuderhfWAdf <- as.data.frame(matrix(data = neuderhfWA, ncol = 1))

neuderlfITdf <- as.data.frame(matrix(data = neuderlfIT, ncol = 1))
neuderlfNCdf <- as.data.frame(matrix(data = neuderlfNC, ncol = 1))
neuderlfWAdf <- as.data.frame(matrix(data = neuderlfWA, ncol = 1))

neuderhfITCH <- str_split_fixed(neuderhfITdf$V1, "_", 2)
neuderhfNCCH <- str_split_fixed(neuderhfNCdf$V1, "_", 2)
neuderhfWACH <- str_split_fixed(neuderhfWAdf$V1, "_", 2)

neuderlfITCH <- str_split_fixed(neuderlfITdf$V1, "_", 2)
neuderlfNCCH <- str_split_fixed(neuderlfNCdf$V1, "_", 2)
neuderlfWACH <- str_split_fixed(neuderlfWAdf$V1, "_", 2)

neuderhfITsigeig2 <- intersect(neuderhfITCH, sigeignam2)
neuderhfNCsigeig2 <- intersect(neuderhfNCCH, sigeignam2)
neuderhfWAsigeig2 <- intersect(neuderhfWACH, sigeignam2)
neuderlfITsigeig2 <- intersect(neuderlfITCH, sigeignam2)
neuderlfNCsigeig2 <- intersect(neuderlfNCCH, sigeignam2)
neuderlfWAsigeig2 <- intersect(neuderlfWACH, sigeignam2)

length(neuderhfITsigeig2)
length(IThighfrqderbadeig2)
length(neuderhfNCsigeig2)
length(NChighfrqderbadeig2)
length(neuderhfWAsigeig2)
length(WAhighfrqderbadeig2)

length(neuderlfITsigeig2)
length(ITlowfrqderbadeig2)
length(neuderlfNCsigeig2)
length(NClowfrqderbadeig2)
length(neuderlfWAsigeig2)
length(WAlowfrqderbadeig2)

citation(package = "limmaDE2")


length(Fstvec)
length(XTXvec)
