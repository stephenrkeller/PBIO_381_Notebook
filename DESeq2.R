# set working directory
setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook")

library("DESeq2")
library("ggplot2")

countsTable <- read.delim("allcountsdataRN_noIT.txt", header=TRUE, stringsAsFactors = TRUE, row.names = 1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_noIT.txt", header=TRUE, stringsAsFactors =  TRUE, row.names = 1)
head(conds)
colData <- as.data.frame(conds)

#########################################

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~devstage + sex + population)
# "population effect" model controlling for differences in devstage and sex

dim(dds)
# [1] 17483    48
dds <- dds[rowSums(counts(dds)) > 1,]

dim(dds)
# [1] 16851    48     Okay, nice, only lost about 600 genes

dds <- DESeq(dds, modelMatrixType = "standard")

resultsNames(dds)
# [1] "Intercept"           "devstage_L3L_vs_AD4" "devstage_PD1_vs_AD4" "devstage_PP1_vs_AD4" "sex_M_vs_F"         
# [6] "population_WA_vs_NC"

res <- results(dds)
str(res)

res <- res[order(res$padj),]
head(res)
# log2 fold change (MLE): population WA vs NC 
# Wald test p-value: population WA vs NC 
# DataFrame with 6 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat       pvalue         padj
# <numeric>      <numeric> <numeric> <numeric>    <numeric>    <numeric>
#   OTAU017482-RA 126.291608     -5.4166340 0.6732174 -8.045892 8.561964e-16 1.268369e-11
# OTAU012716-RA 188.877675      4.2644034 0.5427535  7.856980 3.935069e-15 2.914706e-11
# OTAU008667-RA 231.871115     -0.8736955 0.1267134 -6.895051 5.384538e-12 1.994164e-08
# OTAU012562-RA 251.774364     -0.8774079 0.1270957 -6.903520 5.072966e-12 1.994164e-08
# OTAU013988-RA   4.416955      4.4229857 0.6836393  6.469765 9.815559e-11 2.908154e-07
# OTAU011160-RA  10.241516     -2.5149390 0.4125305 -6.096371 1.085033e-09 2.678947e-06

summary(res)
# out of 16851 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 333, 2% 
# LFC < 0 (down)   : 282, 1.7% 
# outliers [1]     : 85, 0.5% 
# low counts [2]   : 1952, 12% 
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

