#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cloud.r-project.org/")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-d", "--diff"), type="character", default=NULL,
              help="The path to the differential analysis result file (txt) which is needed to restore to homer format", metavar="character"),
  make_option(c("-h", "--homer"), type="character", default=NULL,
              help="The path to the homer quantification file in either rnaQuan or peakQuan folder", metavar="character"),
)
opt_parser = OptionParser("\n\t%prog [options]",
                          option_list=option_list,prog="diff2Homer.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("diff2Homer.R -d .../...vs...txt -h .../rnaQuan/HOMER.rawCount.txt\ndiff2Homer.R -d .../...vs...txt -h .../peakQuan/allRawTags.txt\n", call.=FALSE)
}
opt = parse_args(opt_parser,args)
if(!file.exists(opt$diff)) stop(paste0("The differential result file (",opt$diff,") cannot be loacated!"))
if(!file.exists(opt$homer)) stop(paste0("The homer quantification file (",opt$homer,") cannot be loacated!"))

D <- read.table(opt$diff,sep="\t",header=T,as.is=T,row.names=1,check.names=F,quote="",comment.char="")
H <- read.table(opt$homer,sep="\t",header=T,as.is=T,row.names=1,check.names=F,quote="",comment.char="")

if(sum(rownames(D)%in%rownames(H))<sum(rownames(D)%in%sapply(strsplit(H[,7],"\\|"),head,1))){#RNA
  ix <- sapply(strsplit(H[,7],"\\|"),head,1) %in% rownames(D)
  message("There are ",sum(!ix)," genes not found in differential analysis")
  X <- cbind(H[ix,1:7],D[sapply(strsplit(H[ix,7],"\\|"),head,1),])
}else{#peak
  ix <- rownames(H) %in% rownames(D)
  message("There are ",sum(!ix)," peaks not found in differential analysis")
  X <- cbind(H[ix,1:18],D[rownames(H)[ix],])
}
write.table(X,file=gsub("txt","homer.txt",opt$diff),col.names=NA,quote=F,sep="\t")
