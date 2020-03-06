#!/usr/bin/env Rscript
## peakDiff.R
## initial ----
loadPeakDiff <- function(){
  if(!require(optparse)){
    install.packages("optparse",repos="https://cloud.r-project.org/")
    if(!require(optparse)) stop("Cannot install optparse")
  }
  if(!require(rlang)){
    install.packages("rlang",repos="https://cloud.r-project.org/")
    if(!require(rlang)) stop("Cannot install rlang")
  }
  if(!require(MASS)){
    install.packages("MASS",repos="https://cloud.r-project.org/")
    if(!require(MASS)) stop("Cannot install MASS")
  }
  if(!require(pheatmap)){
    install.packages("pheatmap",repos="https://cloud.r-project.org/")
    if(!require(pheatmap)) stop("Cannot install pheatmap")
  }
  if(!require(DESeq2)){
    BiocManager::install("DESeq2")
    if(!require(DESeq2)) stop("Cannot install DESeq2")
  }
}

suppressWarnings(suppressMessages(loadPeakDiff()))
## input ------------------
args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default].", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="atac", 
              help="The assay type atac, chip and (future: histone) [default= %default], if chip, DESeq library size will be set to sequence depth instead of total tags in peaks", metavar="character"),
  make_option(c("-q", "--qantify"), type="character", #default="", 
              help="The full path to the peak qantification file [default= %default].", metavar="character"),
  make_option(c("-d", "--distal"), type="numeric", default="3000", 
              help="The distance from TSS, peaks outside will be considered [default= %default]", metavar="numeric"),
  make_option(c("-t", "--tss"), type="numeric", #default="3000", 
              help="The distance from TSS, peaks within will be considered. If it is provided, '-d/--distal' will be ignored.", metavar="numeric"),
  make_option(c("-c", "--cutoff"), type="numeric", default="4", 
              help="The cutoff value of tags of a peak [default= %default]", metavar="numeric"),
  make_option(c("-m", "--maxPeak"), type="numeric", default="10000", 
              help="The maximum number of peaks to include in the heatmap [default= %default]", metavar="numeric"),
  make_option(c("-l", "--logFC"), type="numeric", default="1", 
              help="The log fold change to be considered significant [default= %default]", metavar="numeric"),
#  make_option(c("-c", "--homer"), type="character", default="", 
#              help="Additional commands for homer peak annotation (annotatePeaks.pl) [default= %default]", metavar="character"),
  make_option(c("-p", "--padj"), type="numeric", default="0.05", 
              help="The significant p-adj to be considered significant [default= %default]", metavar="numeric")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="peakDiff.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/sample/definition/file is required.\n", call.=FALSE)
}
strSample <- args[1]
opt = parse_args(opt_parser,args[-1])
if(is.null(opt$qantify)){
  print_help(opt_parser)
  stop("The quantification file is required.\n", call.=FALSE)
}
distal <- opt$distal
strOutput <- paste(opt$out,"/",sep="")
if(!dir.exists(strOutput)){
  if(!dir.create(strOutput)){
    stop(paste("Output folder (",strOutput,") cannot be create!\nPlease check the parent folder (",dirname(strOutput),") existence.",sep=""))
  }
}

cat("\n\n=====PLEASE make sure",strSample,"was used to quantify the peaks, or the group information might be mis-matched======\n\n")
## get the sample information----------------
## table 4 (atac) or 5 (chip,histone) columns
## Column 1: grpname
## Column 2: col(rgb)
## Column 3: sampleTag1;sampleTag2;...
## Column 4: sID1;sID2;...
## ChIP & Histome:
## Column 5: (inputTag1;inputTag2;...)
cat("\n\tStep 1: Read the sample group information\n")
tagInfo <- read.table(strSample,sep="\t",row.names=1,as.is = T,comment.char = "")
sID <- COL <- pClass <- c()
for(i in rownames(tagInfo)){
  sID <- c(sID,paste(i,unlist(strsplit(tagInfo[i,3],";")),sep="_"))
  pClass <- c(pClass,rep(i,length(sID)-length(pClass)))
  COL <- c(COL,tagInfo[i,1])
}
names(COL) <- unique(pClass)
print(COL)
if(max(table(pClass))<=1) stop("Replicates are needed!")
## Obtain the distal peaks -----------------
cat("\n\tStep 2: Obtain the tags in peaks\n")
rawTags <- read.table(opt$qantify,as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="",check.names=F)
cat("\t\tTotal peaks:",nrow(rawTags),"\n")
rawC <- as.matrix(rawTags[,-(1:18)])
if(ncol(rawC)!=length(pClass)) stop(paste("ERROR:",opt$qantify,"(",ncol(rawC),") contains different sample number from the sample file(",strSample,")."))
libSize <- as.numeric(sapply(strsplit(sapply(strsplit(colnames(rawC),"\\("),tail,1),"\\."),head,1))
if(opt$assay=="chip"){
  write.table(cbind(idInUse=sID,tagName=basename(sapply(strsplit(colnames(rawC)," "),head,1)),group=pClass,libSize=libSize),
              file=paste(strOutput,"grpInfo.txt",sep=""),sep="\t",row.names=F)
}else{
  write.table(cbind(idInUse=sID,tagName=basename(sapply(strsplit(colnames(rawC)," "),head,1)),group=pClass),file=paste(strOutput,"grpInfo.txt",sep=""),sep="\t",row.names=F)
}
names(libSize) <- colnames(rawC) <- sID

distalC <- rawC[is.na(rawTags$'Distance to TSS')|abs(rawTags$'Distance to TSS')>distal,]
if(!is.null(opt$tss)) distalC <- rawC[!(is.na(rawTags$'Distance to TSS'))|abs(rawTags$'Distance to TSS')<distal,]
distalC <- distalC[apply(distalC,1,function(x){return(sum(x>opt$cutoff))})>1,]
cat("\t\tTotal of",nrow(distalC),"peaks to be considered\n")
## comparison ------------------------------
cat("\n\tStep 3: Use DESeq2 to identify the pair-wised different peak heights\n")
peakDef <- rawTags[rownames(distalC),1:4]
write.table(peakDef,file=paste(strOutput,"/allPeaks.peak",sep=""),
            sep="\t",quote=F,col.names=F)
pheno <- data.frame(row.names = colnames(distalC),grp=pClass)

if(opt$assay=="chip"){
  cat("\t\tsequence depth was set to be the library size normalization for ChIP, \n\t\tthis is consistent with the track\n")
  nf <- matrix(0,nrow=length(libSize),ncol=length(libSize),dimnames = list(names(libSize),names(libSize)))
  diag(nf) <- floor(max(libSize)/1e6)*1e6/libSize
  distalC <- distalC%*%nf
  D <- DESeqDataSetFromMatrix(countData=matrix(as.integer(distalC),nrow=nrow(distalC),dimnames=dimnames(distalC)),
                              colData=pheno,
                              design=as.formula(paste("~",paste(colnames(pheno),collapse="+"))))
  dds <- DESeq(D,betaPrior=TRUE,quiet=T)
  colData(dds)$sizeFactor <- 1
  dds <- DESeq(dds,betaPrior=TRUE,quiet=T)
}else{
  D <- DESeqDataSetFromMatrix(countData=matrix(as.integer(distalC),nrow=nrow(distalC),dimnames=dimnames(distalC)),
                              colData=pheno,
                              design=as.formula(paste("~",paste(colnames(pheno),collapse="+"))))
  dds <- DESeq(D,betaPrior=TRUE,quiet=T)
}

normP <- log2(counts(dds,normalized=T)+1)
## save the results for each pair-wised comparison 
allDBP <- strUpPeak <- c()
strPairwised <- paste(strOutput,"/pairwised/",sep="")
if(!dir.exists(strPairwised)) dir.create(strPairwised)
pdf(paste(strOutput,"/pairwised.pdf",sep=""),width=4,height=4)
par(mar=c(2,2,2,0)+0.2,mgp=c(1,0.1,0),tcl=-0.05)
#imageCOL <- c("#FFFFFFFF","#3300FF","#2D1CFF","#2839FF","#2255FF","#1C71FF","#178EFF","#11AAFF",
#              "#0BC6FF","#06E3FF","#00FFFF","#00FFFF","#17FFE3","#2DFFC6","#44FFAA","#5BFF8E",
#              "#71FF71","#88FF55","#9FFF39","#B5FF1C","#CCFF00",
#             "#CCFF00","#D2F400","#D7E800","#DDDD00","#E3D200","#E8C600","#EEBB00","#F4B000",
#             "#F9A400","#FF9900","#FF9900","#F68800","#EC7700","#E36600","#D95500","#D04400",
#              "#C63300","#BD2200","#B31100","#AA0000")
imageCOL <- c("#FFFFFFFF","#BEBEBE","#B4B4B4","#AAAAAA","#A0A0A0","#969696",
              "#8C8C8C","#828282","#787878","#6E6E6E","#646464","#5A5A5A",
              "#505050","#464646","#3C3C3C","#323232","#282828","#1E1E1E",
              "#141414","#0A0A0A","#000000")
for(i in unique(pClass)){
  peakID <- c()
  cat("\t\tExtract activated peaks for",i,"\n")
  for(j in unique(pClass)){#i:y;j:x
    if(i==j) next
    res <- results(dds,contrast = c("grp",i,j))
    strPair <- paste(strPairwised,j,".vs.",i,".txt",sep="")
    write.table(cbind(data.frame(res),contrast=paste(i,j,sep="-")),
                file=strPair,sep="\t",quote=F,col.names=NA)
    xIndex <- !is.na(res$padj)&res$padj<opt$padj&res$log2FoldChange < -opt$logFC
    yIndex <- !is.na(res$padj)&res$padj<opt$padj&res$log2FoldChange > opt$logFC
    write.table(peakDef[rownames(res)[xIndex],],file=paste(strPair,j,".vs.",i,"_",j,".peak",sep=""),sep="\t",quote=F,col.names=NA)
    write.table(peakDef[rownames(res)[yIndex],],file=paste(strPair,j,".vs.",i,"_",i,".peak",sep=""),sep="\t",quote=F,col.names=NA)
    
    Col <- rep("gray",nrow(res))
    Col[xIndex] <- COL[j]
    Col[yIndex] <- COL[i]
    x <- apply(normP[rownames(res),pClass==j,drop=F],1,mean)
    y <- apply(normP[rownames(res),pClass==i,drop=F],1,mean)
    ylim <- xlim <- range(c(x,y))
    plot(c(),c(),xlab=j,ylab=i,xlim=xlim,ylim=ylim,main="log2 Mean normalized tag")
    for(k in c("gray",COL[c(i,j)])){
      index <- rep(T,sum(Col==k))
      if(k == "gray"){
        tryM <- try(f1 <- MASS::kde2d(x[Col==k],y[Col==k],n=100),silent=T)
        if(!is.null(names(tryM))){
          image(f1,col=imageCOL,add=T)
          imageZero <- diff(range(f1$z))/length(imageCOL)
          index <- apply(cbind(x[Col==k],y[Col==k]),1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
        }
      }
      points(x[Col==k][index],y[Col==k][index],pch=20,col=k,cex=1,useDingbats = F)
    }
    lines(range(c(xlim,ylim)),range(c(xlim,ylim)),col="gray")
    lines(range(c(xlim,ylim)),range(c(xlim,ylim))+opt$logFC,col="gray",lty=2)
    lines(range(c(xlim,ylim)),range(c(xlim,ylim))-opt$logFC,col="gray",lty=2)
    mtext(paste(i,": ",sum(yIndex),sep=""),3,line=-1,adj=0.1,col=COL[i])
    mtext(paste(j,": ",sum(xIndex),sep=""),4,line=-1,adj=0.1,col=COL[j])
    
    peakID <- c(peakID,rownames(res)[yIndex])
  }
  peakID <- unique(peakID)
  if(length(peakID)>100) write.table(peakDef[peakID,],file=paste(strOutput,"/",i,".act.peak",sep=""),
                                     sep="\t",quote=F,col.names=F)
  allDBP <- c(allDBP,peakID)
}
a <- dev.off()
allDBP <- unique(allDBP)
if(length(allDBP)>10){
  write.table(peakDef[allDBP,],file=paste(strOutput,"/all.act.peak",sep=""),
              sep="\t",quote=F,col.names=F)
  if(length(allDBP)>opt$maxPeak) cat("\t\tSubset (",opt$maxPeak,") of total",length(allDBP),"peaks were selected for heatmap\n")
  subDBP <- allDBP[sample(length(allDBP),min(opt$maxPeak,length(allDBP)))]
  pdf(paste(strOutput,"/all.act.pdf",sep=""),height = 9)#,onefile=FALSE
  heatCol <- c("#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
               "#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5","#F8EFEA","#FBE0D1",
               "#FAC9B1","#F5AD8C","#E88B6E","#D96752","#C6413E","#B31B2C",
               "#8E0C25","#67001F")
  pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),labels_row=rep("",length(subDBP)),
           color = heatCol,
           annotation_col = data.frame(row.names=colnames(normP),
                                       grp=pClass))
  pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),scale="row",labels_row=rep("",length(subDBP)),
           color = heatCol,
           annotation_col = data.frame(row.names=colnames(normP),
                                       grp=pClass))
  pheatmap(normP[subDBP,],annotation_colors=list(grp=COL),scale="row",labels_row=rep("",length(subDBP)),
           color = heatCol,cluster_cols=F,
           annotation_col = data.frame(row.names=colnames(normP),
                                       grp=pClass))
  a <- dev.off()
}

## find motifs for activated peaks ----------------------------
cat("\tStep 4: Find motifs for activated peaks\n")
for(i in list.files(strOutput,"act.peak",full.names=T)){
  if(grepl("all\\.act\\.peak",i)) next
  cat("\t\tFind motifs for",basename(i),"\n")
  strCMD <- paste("findMotifsGenome.pl",i,opt$genome,gsub("peak$","bgR",i))
  cat("\t\t",strCMD,"\n")
  system(paste(strCMD,"2>/dev/null"))
}
## plot top 5 of motifs ----------------------
cat("\tStep 5: plot top 5 motifs \n")
conn <- file(paste(strOutput,"/motif",sep=""),"w")
cat("motif\tcol\tknown\tdenovo\t\n",file=conn)
for(i in list.dirs(strOutput,recursive=F)){
  if(!grepl("bgR$",i)) next
  cat(i,COL[gsub("\\.act\\.bgR","",basename(i))],
      paste(1:5,collapse=","),
      paste(1:5,collapse=","),
      "\n",sep="\t",file=conn)
}
close(conn)
strCMD <- paste("plotMotif.R",paste(strOutput,"/motif",sep=""))
cat("\t\t",strCMD,"\n")
system(paste(strCMD,"2>/dev/null"))

## ---------------------
cat("\nDifferential analysis is done successfully!\n")
















