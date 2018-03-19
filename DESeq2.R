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

res_pop <- results(dds, name="population_WA_vs_NC", alpha=0.05)

res_pop <- res_pop[order(res_pop$padj),]
summary(res_pop)

######################## Data visualization
plotMA(res_pop, main="DESeq2",ylim=c(-2,2))
abline(h=c(-1,1), col="blue", lwd=2)

# sex effect?

res_sex <- results(dds, name="sex_M_vs_F", alpha=0.5)
plotMA(res_sex, main="DESeq2",ylim=c(-2,2))
abline(h=c(-1,1), col="blue", lwd=2)

######## PCA

vsd <- vst(dds,blind=FALSE)

data <- plotPCA(vsd,intgroup=c("population","devstage","sex"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
data$devstage <- factor(data$devstage, levels=c("L3L","PP1","PD1","AD4"), labels = c("L3L","PP1","PD1","AD4"))

pdf(file="PCA_sex_by_stage.pdf",height=5.5,width=5.5)
ggplot(data,aes(PC1,PC2,color=sex,shape=devstage)) +
  geom_point(size=4,alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()
dev.off()

ggplot(data,aes(PC1,PC2,color=population,shape=devstage)) +
  geom_point(size=4,alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

############ the next several lines will make a cluster heatmap of all the samples vs all samples

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$population, vsd$devstage, vsd$sex, sep="-")
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette(rev(brewer.pal(9,"Purples")))(255) # 255 delimits the number of bins in the scale

colors <- colorRampPalette(c("Purple", "Green"))(255)

library("pheatmap")

pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists, col=colors)

############### Let's look at individual genes!!

d <- plotCounts(dds, gene="OTAU017482-RA", intgroup= (c("population","sex","devstage")), returnData=TRUE)
d

p <- ggplot(d, aes(x=devstage, y=count, shape=sex, colour = population)) + 
  theme_minimal() + theme(text=element_text(size=20), panel.grid.major = element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) + 
  scale_x_discrete(limits=c("L3L","PP1","PD1","AD4"))
p

############### Save outputs

write.csv(res_pop, file="DGE_NCvsWA_pop_classDay6.csv", row.names=TRUE, quote=FALSE)

########## change p-value to -log(pvalue) and export as .csv with row names

neglogpval <- as.matrix(-log(res_pop$pvalue))
head(neglogpval)
-log(8.561964e-16) # checking first value

res_pop_negpval <- cbind(row.names(res_pop),neglogpval)
head(res_pop_negpval)

colnames(res_pop_negpval)=c("gene","neglogpval")

write.csv(res_pop_negpval, file="DGE_NCvsWA_pop_neglogpval_classDay6.csv", row.names=FALSE,quote=FALSE,col.names=TRUE)

# ----------------------------- Day 7 Transcriptomics

colData$group <- factor(paste0(colData$population, "-", colData$devstage, "-", colData$sex))
head(colData)

dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design= ~ group)
dds <- dds[rowSums(counts(dds)) > 1, ]
dim(dds)

dds <- DESeq(dds, modelMatrixType = "standard", parallel=T)

resultsNames(dds)

# [1] "Intercept"                  "group_NC.AD4.M_vs_NC.AD4.F" "group_NC.L3L.F_vs_NC.AD4.F" "group_NC.L3L.M_vs_NC.AD4.F"
# [5] "group_NC.PD1.F_vs_NC.AD4.F" "group_NC.PD1.M_vs_NC.AD4.F" "group_NC.PP1.F_vs_NC.AD4.F" "group_NC.PP1.M_vs_NC.AD4.F"
# [9] "group_WA.AD4.F_vs_NC.AD4.F" "group_WA.AD4.M_vs_NC.AD4.F" "group_WA.L3L.F_vs_NC.AD4.F" "group_WA.L3L.M_vs_NC.AD4.F"
# [13] "group_WA.PD1.F_vs_NC.AD4.F" "group_WA.PD1.M_vs_NC.AD4.F" "group_WA.PP1.F_vs_NC.AD4.F" "group_WA.PP1.M_vs_NC.AD4.F"

res_pop_PP1_F <- results(dds, contrast = list(c("group_NC.PP1.F_vs_NC.AD4.F"),c("group_WA.PP1.F_vs_NC.AD4.F")), listValues=c(1/2,-1/2),alpha=0.05)

res_pop_PP1_F <- res_pop_PP1_F[order(res_pop_PP1_F$padj),]
head(res_pop_PP1_F)

summary(res_pop_PP1_F)

# out of 16851 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)     : 46, 0.27% 
# LFC < 0 (down)   : 31, 0.18% 
# outliers [1]     : 17, 0.1% 
# low counts [2]   : 2940, 17% 
# (mean count < 3)

######## Pull out significant genes to be able to make a heatmap

sig_pop_PP1_F <- res_pop_PP1_F[which(res_pop_PP1_F$padj < 0.05), ] # great way to subset data
dim(sig_pop_PP1_F) # check

write.csv(sig_pop_PP1_F, file="Gene_list")

sig_pop_PP1_F_df <- as.data.frame(sig_pop_PP1_F)
sig_pop_PP1_F_df$Row.names <- rownames(sig_pop_PP1_F_df)
dim(sig_pop_PP1_F_df)

genesOfInterest_pop_PP1_F <- c(sig_pop_PP1_F_df$Row.names)
length(genesOfInterest_pop_PP1_F)

vsd <- vst(dds, blind=FALSE)

dds$combined = factor(paste0(dds$population, "-", dds$devstage, "-", dds$sex))
dds$combined <- factor(dds$combined, levels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"), labels=c("WA-L3L-F","WA-L3L-M","WA-PP1-F","WA-PP1-M","WA-PD1-F","WA-PD1-M","WA-AD4-F","WA-AD4-M","NC-L3L-F","NC-L3L-M","NC-PP1-F","NC-PP1-M","NC-PD1-F","NC-PD1-M","NC-AD4-F","NC-AD4-M"))

baseMeanPerGrp <- sapply( levels(dds$combined), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$combined == lvl] ) )

head(baseMeanPerGrp)
dim(baseMeanPerGrp)

# pulls out normalized counts (avg of 3 reps) for all of our significant genes
m <- baseMeanPerGrp[genesOfInterest_pop_PP1_F, c("WA-PP1-F", "WA-PP1-M", "NC-PP1-F", "NC-PP1-M")]

head(m)
dim(m)

mat_scaled = t(apply(m, 1, scale))
head(mat_scaled)

library(pheatmap)

pheatmap(mat_scaled, labels_col=c("WA-PP1-F", "WA-PP1-M", "NC-PP1-F", "NC-PP1-M"), cluster_cols=FALSE, cluster_rows=TRUE)

#######################

# Export counts to use as the input for WGCNA

norm.counts <- counts(dds, normalized=TRUE)
dim(norm.counts)

write.csv(norm.counts, file="beetle_norm_counts.csv", row.names=TRUE, quote=FALSE)
