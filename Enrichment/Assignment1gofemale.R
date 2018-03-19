
# GO_MWU uses continuous measure of significance (such as fold-change or -log(p-value) ) to identify GO categories that are significantly enriches with either up- or down-regulated genes. The advantage - no need to impose arbitrary significance cutoff.

# If the measure is binary (0 or 1) the script will perform a typical "GO enrichment" analysis based Fisher's exact test: it will show GO categories over-represented among the genes that have 1 as their measure. 

# On the plot, different fonts are used to indicate significance and color indicates enrichment with either up (red) or down (blue) regulated genes. No colors are shown for binary measure analysis.

# The tree on the plot is hierarchical clustering of GO categories based on shared genes. Categories with no branch length between them are subsets of each other.

# The fraction next to GO category name indicates the fracton of "good" genes in it; "good" genes being the ones exceeding the arbitrary absValue cutoff (option in gomwuPlot). For Fisher's based test, specify absValue=0.5. This value does not affect statistics and is used for plotting only.

# Stretch the plot manually to match tree to text

# Mikhail V. Matz, UT Austin, February 2015; matz@utexas.edu

################################################################
# First, press command-D on mac or ctrl-shift-H in Rstudio and navigate to the directory containing scripts and input files. Then edit, mark and execute the following bits of code, one after another.
setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook/Enrichment")

# Can run for all 6 inputs, mean and max for each of the three pop. pair FSTs and can do for each GO division BP, MF, CC

# Edit these to match your data file names: 
input="DGE_female_neglogpval.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
## Make sure this is saved as Unix (LF) - open in TextWrangler, save as, change from Classic Mac to Unix (LF)!!!
goAnnotations="gene_annotation_only.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download from http://www.geneontology.org/GO.downloads.ontology.shtml
goDivision="BP" # either MF, or BP, or CC
source("gomwu.functions.R")


# Calculating stats. It takes ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. 
           Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# go.obo gene_annotation_only.tab DGE_NCvsWA_pop_pval.csv BP largest=0.1 smallest=5cutHeight=0.25

# Run parameters:

# largest GO category as fraction of all genes (largest)  : 0.1
# smallest GO category as # of genes (smallest)  : 5
# clustering threshold (clusterCutHeight) : 0.25

# -----------------
# retrieving GO hierarchy, reformatting data...

# -------------
# go_reformat:
# Genes with GO annotations, but not listed in measure table: 5

# Terms without defined level (old ontology?..): 0
# -------------
# -------------
# go_nrify:
# 1568 categories, 1294 genes; size range 5-129.4
# 34 too broad
# 1131 too small
# 403 remaining

# removing redundancy:

# calculating GO term similarities based on shared genes...
# 209 non-redundant GO categories of good size
# -------------

# Secondary clustering:
# calculating similarities....
# Continuous measure of interest: will perform MWU test
# 3  GO terms at 10% FDR  for DGE_NCvsWA_pop_pval - 3 very different, interesting, transport, signal transduction, and metabolism

# Continuous measure of interest: will perform MWU test
# 3  GO terms at 10% FDR  for DGE_NCvsWA_pop_absstat -  same 3 types of categories as above
# GO terms dispayed:  3 
# "Good genes" accounted for:  61 out of 393 ( 16% )  - but here this is different, above showed 0 "good" genes, so somthing about pvals vs the actual statistic

# Continuous measure of interest: will perform MWU test
# 5  GO terms at 10% FDR  for DGE_NCvsWA_pop_absLFC
# GO terms dispayed:  5 - all related to metabolism - LFC is not as good of metric to use.
# "Good genes" accounted for:  0 out of 0 ( NaN% )

###################### Plotting results
library(ape)
quartz()
gomwuPlot(input,goAnnotations,goDivision,
          absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.5 if you are doing Fisher's exact test for standard GO enrichment.
          level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
          level2=0.05, # FDR cutoff to print in regular (not italic) font.
          level3=0.01, # FDR cutoff to print in large bold font.
          txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
          treeHeight=0.5, # height of the hierarchical clustering tree
          #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
)
# manually rescale the plot so the tree matches the text 
# if there are too many categories displayed, try make it more stringent with level1=0.01,level2=0.005,level3=0.001.  