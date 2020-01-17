#!/usr/bin/env Rscript
############################
## plotMergeVenn.R
##
############################
if(!require(Vennerable,quiet=T)){
  if(!require(RBGL,quiet=T)) BiocManager::install("RBGL")
  if(!require(graph,quiet=T)) BiocManager::install("graph")
  devtools::install_github("js229/Vennerable")
}

if(!require(Vennerable,quiet=T)) stop("Cannot install package 'Vennerable'!")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("usage: plotMergeVenn.R /path/to/homer/mergePeak/venn/file\n")
  cat("e.g.: plotMergeVenn.R mergedVenn\n")
  q()
}
strF <- args[1]
#strF <- "~/scratch/YA/H3k27ac_NCorWT/DCA/withATAC/ATAC_H3K27ac.venn"#args[1]#
if(file.exists(strF)){
  A <- read.table(strF,sep="\t",header=T,check.names=F,as.is=T)
  colnames(A) <- gsub("_idr.txt","",basename(colnames(A)))
  setNames <- colnames(A)[1:(grep("^Total$",colnames(A))-1)]
  ids <- w <- c()
  for(i in 1:nrow(A)){
    tmp <- rep(0,length(setNames))
    tmp[nchar(A[i,1:length(setNames)])>0] <- 1
    ids <- c(ids,paste(tmp,collapse=""))
    w <- c(w,A[i,"Total"])
  }
  names(w) <- ids
  V3a <- Venn(SetNames=setNames,Weight=w)
  Weights(V3a)[is.na(Weights(V3a))] <- 0
  #print(V3a)
  pdf(paste(gsub("\\.txt","",strF),".pdf",sep=""))
  if(length(setNames)>2)
    plot(V3a, type = "ChowRuskey")
  else
    plot(V3a)
  tmp <- dev.off()
}