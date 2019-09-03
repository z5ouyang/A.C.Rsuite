#!/usr/bin/env Rscript
############################
## plotMergeVenn.R
##
############################
if(!require(Vennerable,quiet=T)) install.packages("Vennerable")
if(!require(Vennerable,quiet=T)) stop("Cannot install package 'Vennerable'!")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("usage: plotMergeVenn.R /path/to/homer/mergePeak/venn/file\n")
  cat("e.g.: plotMergeVenn.R mergedVenn\n")
  q()
}
strF <- args[1]#"~/scratch/Lung/Results/IP/ATAC/IDR/mergeIDRs/IMmono.txt"#
if(file.exists(strF)){
  A <- read.table(strF,sep="\t",header=T,check.names=F,as.is=T)
  colnames(A) <- gsub("_idr.txt","",basename(colnames(A)))
  setNames <- colnames(A)[1:(grep("^Total$",colnames(A))-1)]
  ids <- w <- c()
  for(i in 1:nrow(A)){
    tmp <- rep(0,length(setNames))
    tmp[nchar(A[i,1:length(setNames)])>0] <- 1
    ids <- c(ids,paste(tmp,collapse=""))
    w <- c(w,A[,"Total"])
  }
  names(w) <- ids
  V3a <- Venn(SetNames=setNames,Weight=w)
  #print(V3a)
  pdf(paste(gsub("\\.txt","",strF),".pdf",sep=""))
  if(length(setNames)>2)
    plot(V3a, type = "ChowRuskey")
  else
    plot(V3a)
  tmp <- dev.off()
}