setwd("~/Documents/UVM_2018/PBIO381/PBIO_381_Notebook")

install.packages(c("Cairo","ggplot2","gridExtra","gtable","tidyr","devtools"),dependencies=T)
devtools::install_github('royfrancis/pophelper')

library(pophelper)


admixfiles=list.files(path=("ADMIXTURE/"),full.names=T)
admixlist=readQ(files=admixfiles,filetype="basic")

# metadata: sample id and pop from beetle.pop file
metadata=read.table("cols_data.txt",header=T)

# format metadata to a data frame and ind variables as chars. for plotting
metadata2=data.frame(sampleid=metadata[,1], population=metadata[,2])

metadata2$sampleid=as.character(metadata2$sampleid)
metadata2$population=as.character(metadata2$population)

# add in the sample id to the different Q files for plotting
if(length(unique(sapply(admixlist,nrow)))==1)
  admixlist <- lapply(admixlist,"rownames<-",metadata2$sampleid)

# three is the k = 4
head(admixlist[[3]])

p <- plotQ(admixlist[c(1,2,3)],
           returnplot=T,exportplot=T,quiet=T,basesize=11, imgoutput="join", 
           showyaxis=T, showticks=T, panelspacer=0.4, useindlab=F, showindlab=F, 
           grplab=metadata2[2], grplabsize=4, linesize=1, barsize=1, pointsize=3, 
           panelratio=c(4,1), divsize = 0.75, pointcol="white", showtitle=T, 
           titlelab="ADMIXTURE analysis on O. tauri, SNPs thinned to 1 per kb", 
           splab=c("K=2","K=3","K=4"), outputfilename="ADMIXTURE_Otauri",
           imgtype="pdf",height=3,width=25)

plot(p$plot[[1]])
