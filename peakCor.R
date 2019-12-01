#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(MASS)))) install.packages("MASS",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(pheatmap)))) install.packages("pheatmap",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse) || !require(MASS) || !require(pheatmap))
  stop("R packages of pheatmap, MASS or optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="atac", 
              help="The assay type atac, chip and (future: histone) [default= %default]", metavar="character"),
  make_option(c("-c", "--homer"), type="character", default=" -minDist 200 -size 200", 
              help="Additional commands for homer peak calling [default= '%default']", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="peakCor.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/sample/definition/file is required.\n", call.=FALSE)
}
strSample <- args[1]
opt = parse_args(opt_parser,args[-1])
strOutput <- paste(opt$out,"/",sep="")
if(!dir.exists(strOutput)){
  if(!dir.create(strOutput)){
    stop(paste("Output folder (",strOutput,") cannot be create!\nPlease check the parent folder (",dirname(strOutput),") existence.",sep=""))
  }
}
opt$homer <- trimws(opt$homer)
## processing ------
## table 4 (atac) or 5 (chip,histone) columns
## Column 1: grpname
## Column 2: col(rgb)
## Column 3: sampleTag1;sampleTag2;...
## Column 4: sID1;sID2;...
## ChIP & Histome:
## Column 5: (inputTag1;inputTag2;...)
exps <- read.table(strSample,sep="\t",comment.char="",as.is=T)
strTmp <- paste(strOutput,"peakCor_tmp",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)
limFun <- function(x){
  step <- diff(range(x))/20
  return(range(x)-c(step,-step))
}
heatCol <- c("#99000D","#A80712","#B80F17","#C8161C","#D42120","#DF2C25",
             "#EB372A","#F14432","#F5533B","#F96245","#FB7050","#FB7C5C",
             "#FB8969","#FC9676","#FCA385","#FCB094","#FCBDA3","#FCCAB5",
             "#FDD7C7","#FEE5D9")
imageCOL <- c("#FFFFFFFF","#3300FF","#2D1CFF","#2839FF","#2255FF","#1C71FF","#178EFF","#11AAFF",
              "#0BC6FF","#06E3FF","#00FFFF","#00FFFF","#17FFE3","#2DFFC6","#44FFAA","#5BFF8E",
              "#71FF71","#88FF55","#9FFF39","#B5FF1C","#CCFF00",
              "#CCFF00","#D2F400","#D7E800","#DDDD00","#E3D200","#E8C600","#EEBB00","#F4B000",
              "#F9A400","#FF9900","#FF9900","#F68800","#EC7700","#E36600","#D95500","#D04400",
              "#C63300","#BD2200","#B31100","#AA0000")
peakSumm <- c()
for(i in 1:nrow(exps)){
  grpID <- exps[i,1]
  cat("\tAnalysing ",grpID,"\n",sep="")
  tags <- unlist(strsplit(exps[i,3],";"))
  sIDs <- paste(grpID,unlist(strsplit(exps[i,4],";")),sep="_")
  if(opt$assay=="chip"){
    inputs <- unlist(strsplit(exps[i,5],";"))
  }
  if(length(tags)!=length(sIDs) || (exists("inputs")&&length(tags)!=length(inputs)))
    stop("the number of samples is different from the number of sample id or inputs!")
  strPeak <- c()
  for(j in 1:length(tags)){
    if(!dir.exists(tags[j])){
      stop(paste(tags[j],"does NOT exist!"))
    }
    strPeak <- c(strPeak,paste(strTmp,"/",sIDs[j],".peaks",sep=""))
    if(opt$assay=="atac"){
      strCMD <- paste("findPeaks",tags[j],"-style factor",opt$homer,"-o",tail(strPeak,1))
    }else{
      strCMD <- paste("findPeaks",tags[j],"-style factor -i",inputs[j],opt$homer,"-o",tail(strPeak,1))
    }
    cat("\n\tStep 1: Finding peaks for a sample",sIDs[j],"\n")
    cat("\t\t",strCMD,"\n")
    #if(!file.exists(tail(strPeak,1)))
      system(paste(strCMD,"2>/dev/null"))
    ## reading the peak information
    key <- "sID"
    value <- sIDs[j]
    
    for(k in scan(tail(strPeak,1),character(),sep="\n",nlines=30,quiet=T)){
      if(grepl("=",k)){
        res <- unlist(strsplit(gsub("#","",k),"="))
        key <- c(key,trimws(res[1]))
        value <- c(value,trimws(res[2]))
      }
    }
    names(value) <- key
    peakSumm <- rbind(peakSumm,value)
  }
  if(length(tags)<2) next
  
  cat("\n\tStep 2: Mergeing all peaks from replicates for",grpID,"\n")
  strMerge <- paste(strTmp,"/",grpID,".peaks",sep="")
  strCMD <- paste("mergePeaks",paste(strPeak,collapse = " "),">",strMerge)
  cat("\t\t",strCMD,"\n")
  #if(!file.exists(strMerge))
    system(paste(strCMD,"2>/dev/null"))

  cat("\n\tStep 3: Obtaining the normalized tag counts of all replicates on the merged replicate peaks for",grpID,"\n")
  strAnno <- gsub("peaks$","anno",strMerge)
  strCMD <- paste("annotatePeaks.pl",strMerge,opt$genome,"-d",paste(tags,collapse=" "),">",strAnno)
  cat("\t\t",strCMD,"\n")
  #if(!file.exists(strAnno))
    system(paste(strCMD,"2>/dev/null"))

  cat("\n\tStep 4: Reading the tag counts and ploting the correlation among replicates of",grpID,"\n")
  X <- read.table(strAnno,sep="\t",header=T,as.is=T,row.names=1,comment.char = "",quote="",check.names=F)
  reads <- log2(X[,19:ncol(X)]+1)
  colnames(reads) <- sIDs
  strCor <- paste(strOutput,grpID,"_cor.pdf",sep="")
  pdf(strCor,width=4,height=4)
  par(cex=0.5,mgp=c(1.5,0.5,0),mar=c(2.5,2.5,0,0)+0.5)
  COR <- matrix(1,nrow=ncol(reads),ncol=ncol(reads),dimnames=list(sIDs,sIDs))
  for(j in 1:(ncol(reads)-1)){
    for(k in (j+1):ncol(reads)){
      index <- apply(reads[,c(j,k)],1,max)>0
      tmp <- reads[index,c(j,k)]
      tryM <- try(f1 <- kde2d(tmp[,1],tmp[,2],n=500),silent=T)
      if(is.null(names(tryM))){
        plot(c(),c(),xlab=sIDs[j],ylab=sIDs[k],xlim=limFun(f1$x),ylim=limFun(f1$y))
        index <- rep(T,nrow(tmp))
      }else{
        image(f1,col=imageCOL,xlab=sIDs[j],ylab=sIDs[k],xlim=limFun(f1$x),ylim=limFun(f1$y))
        imageZero <- diff(range(f1$z))/length(imageCOL)
        index <- apply(tmp,1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
      }
      points(tmp[index,1],tmp[index,2],col=imageCOL[2],pch=20)
      lines(range(c(limFun(f1$x),limFun(f1$y))),range(c(limFun(f1$x),limFun(f1$y))),col="gray")
      COR[j,k] <- COR[k,j] <- cor(tmp[,1],tmp[,2])
      legend("topleft",paste("r=",round(COR[j,k],3)),bty="n")
    }
  }
  print(COR)
  pheatmap(COR,rev(heatCol),cellwidth=20,cellheight=20,
           cluster_rows=F,cluster_cols=F,
           fontsize=8,
           display_numbers=T,number_color="black")
  a <- dev.off()
}
cat("\n\tStep 5: Obtaining all each individual peak information\n")
write.table(peakSumm,file=paste(strOutput,"peak_summary.txt",sep=""),sep="\t",row.names=F,quote=F)


cat("\nThe correlation process is completed successfully!\n")