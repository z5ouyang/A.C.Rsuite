#!/usr/bin/env Rscript
############################################
## alignStats.R
## 
#############################################
if(!require(optparse,quietly = T)) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse))
  stop("R packages of optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
#  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
#              help="Output directory path [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file > path/tp/result/folder/stat.txt",
                          option_list=option_list,prog="alignStats.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/sample/definition/file is required.\n", call.=FALSE)
}
#cat("test\n")
strSample <- args[1]
if(dir.exists(strSample)){
  strDir <- list.dirs(strSample,recursive = F)
  names(strDir) <- basename(strDir)
}else{
  strDir <- c()
  exps <- read.table(strSample,sep="\t",comment.char="",as.is=T)
  for(i in 1:nrow(exps)){
    grpID <- exps[i,1]
    one <- unlist(strsplit(exps[i,3],";"))
    oneN <- paste(grpID,unlist(strsplit(exps[i,4],";")),sep="_")
    if(length(one)!=length(oneN)){
      stop(paste("Number of Tag directories (",length(one),") is different from sample ID name (",length(oneN),") for ",grpID,sep=""))
    }
    names(one) <- oneN
    strDir <- c(strDir,one)
  }
  
}

#cat("Extract alignment stats from:\n",paste(strDir,collapse="\n"),"\n",sep="\t")
######################################################
## stats
stat.names <- c("alignerTotal","alignerUnique","alignerUniqueRate","alignerMulti","homerUniPos","homerTotal","tagPosition","FragLength","peakSize","homerAvgLength","mitoNum","mitoRate")
stat <- matrix(NA,nrow=length(strDir),ncol=length(stat.names),dimnames=list(names(strDir),stat.names))

for(i in names(strDir)){
  #cat(i,"\n")
  strLog <- list.files(strDir[i],"star.log$",full.names=T)
  if(length(strLog)==1){
    sStat <- scan(strLog,character(),sep="\t",quiet=T)
    sStat <- as.numeric(sStat[c(grep("input reads",sStat)+1,grep("Uniquely mapped reads number",sStat)+1,grep("Number of reads mapped to multiple loci",sStat)+1)])
    stat[i,1:4] <- c(sStat[1:2],sStat[2]/sStat[1],sStat[3])
  }
  strLog <- list.files(strDir[i],"bowtie2.log$",full.names=T)
  if(length(strLog)==1){
    sStat <- scan(strLog,character(),sep=" ",quiet=T)
    sStat <- as.numeric(sStat[c(grep("reads;",sStat)-1,
                                grep("exactly",sStat)-3,
                                grep(">1",sStat)-3)])
    stat[i,1:4] <- c(sStat[1:2],sStat[2]/sStat[1],sStat[3])
  }
  res <- read.table(paste(strDir[i],"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T)
  stat[i,"homerUniPos"] <- as.numeric(res[1,"Unique.Positions"])#grepl("genome=",res[,1])
  stat[i,"homerTotal"] <- as.numeric(res[1,"Total.Tags"])#grepl("genome=",res[,1])
  stat[i,"tagPosition"] <- res[1,"Total.Tags"]/res[1,"Unique.Positions"]
  stat[i,"FragLength"] <- as.numeric(gsub("fragmentLengthEstimate=","",res[grepl("fragmentLengthEstimate",res[,1]),1]))
  stat[i,"peakSize"] <- as.numeric(gsub("peakSizeEstimate=","",res[grepl("peakSizeEstimate",res[,1]),1]))
  stat[i,"homerAvgLength"] <- as.numeric(gsub("averageTagLength=","",res[grepl("averageTagLength",res[,1]),1]))
  if(file.exists(paste(strDir[i],"/tagInfo_with_M.txt",sep=""))) res <- read.table(paste(strDir[i],"/tagInfo_with_M.txt",sep=""),sep="\t",as.is=T,header=T)
  index <- res[,1]=="chrM"
  if(sum(index)==1){
    stat[i,"mitoNum"] <- as.numeric(res[index,"Total.Tags"])
    stat[i,"mitoRate"] <- stat[i,"mitoNum"]/ as.numeric(res[1,"Total.Tags"])
  }else if(sum(index)==0){
    stat[i,"mitoRate"] <- stat[i,"mitoNum"] <- 0
  }else{
    stop("ERROR: more than one chrM")
  }
}
cat("\t",paste(colnames(stat),collapse="\t"),"\n",sep="")
for(i in rownames(stat)) cat(i,"\t",paste(stat[i,],collapse="\t"),"\n",sep="")

#write.table(stat,file=paste(strOutput,"/HOMER.stats.txt",sep=""),sep="\t",quote=F,col.names=NA)
