#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(rlang)))) install.packages("rlang",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(DESeq2)))) BiocManager::install("DESeq2")
if(!suppressWarnings(suppressMessages(require(pheatmap)))) install.packages("pheatmap",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(MASS)))) install.packages("MASS",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(plotrix)))) install.packages("plotrix",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse) || !require(DESeq2) || !require(pheatmap) || !require(MASS) || !require(plotrix))
  stop("R packages of optparse, pheatmap, MASS, plotrix or DESeq2 cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-c", "--counts"), type="character", #default="mm10", 
              help="required. The path to the raw counts file which only contains sample columns.", metavar="character"),
  make_option(c("-t", "--tpm"), type="character", #default="250", 
              help="required. The path to the TPM file which only contains sample columns.", metavar="character"),
  make_option(c("-m", "--minTPM"), type="numeric", default="8", 
              help="The minumal TPM required in at least 2 samples", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="rnaDiff.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/sample/definition/file is required.\n", call.=FALSE)
}
strSample <- args[1]
opt <- parse_args(opt_parser,args[-1])
if (is.null(opt$counts)||is.null(opt$tpm)){
  print_help(opt_parser)
  stop("Counts and TPM are both required.\n", call.=FALSE)
}
strOutput <- paste(opt$out,"/",sep="")
if(!dir.exists(strOutput)){
  if(!dir.create(strOutput)){
    stop(paste("Output folder (",strOutput,") cannot be create!\nPlease check the parent folder (",dirname(strOutput),") existence.",sep=""))
  }
}
## processing ------
## table 4 columns
## Column 1: grpname
## Column 2: col(rgb)
## Column 3: sampleTag1;sampleTag2;...
## Column 4: sID1;sID2;...
## obtain the sample group information ------
cat("\n\tStep 1: check the group information\n")
sInfo <- read.table(strSample,as.is=T,comment.char="",sep="\t",row.names = 1)#
COL <- sInfo[,1]
names(COL) <- rownames(sInfo)
sID <- pClass <- c()
for(i in rownames(sInfo)){
  one <- paste(i,unlist(strsplit(sInfo[i,3],";")),sep="_")
  pClass <- c(pClass,rep(i,length(one)))
  sID <- c(sID,one)
}
rawT <- as.matrix(read.table(opt$tpm,sep="\t",header=T,row.names=1))
rawC <- as.matrix(read.table(opt$counts,sep="\t",header=T,row.names=1))
grp <- cbind(TPMid=colnames(rawT),
             CountID=colnames(rawC),
             sID=sID,
             grp=pClass)
print(grp)
## pair-wised comparisons -----------------
cat("\n\tStep 2: pair-wised group comparisons\n")
DESeq2.scatter <- function(X,pClass,strPath,COL,cutOff=4,showX=NULL,logFC=1,padj=0.05,
                           title="",pSize=1,prefix=""){
  com <- names(COL)
  if(length(com)!=2) stop("The specified category number is NOT 2.")
  if(sum(table(pClass)<2)){
    cat("\Ignore comparison of",paste(com,collapse = ".vs."),"for less than 2 replicates\n")
    return(c())
  }
  if(!is.null(showX)) Data <- X[apply(showX,1,function(x){sum(x>cutOff)})>1,]
  else Data <- X[apply(X,1,function(x){sum(x>cutOff)})>1,]
  
  pheno <- data.frame(row.names=colnames(Data),grp=pClass)
  D <- DESeqDataSetFromMatrix(countData=matrix(as.integer(round(Data)),nrow=nrow(Data),dimnames=dimnames(Data)),
                              colData=pheno,
                              design=as.formula(paste("~",paste(colnames(pheno),collapse="+"))))
  dds <- DESeq(D,betaPrior=TRUE,quiet=T)
  if(!is.null(showX)) Data <- log2(showX[rownames(Data),colnames(Data)]+1)
  else Data <- log2(counts(dds,normalized=T)+1)
  if(nchar(prefix)>0) prefix <- paste(prefix,"_",sep="")
  strF <- paste(strPath,"/",prefix,paste(rev(com),collapse=".vs."),".scatter.pdf",sep="")
  
  ## save the DESeq2 diff results
  res <- cbind(data.frame(results(dds,contrast=c("grp",com))),contrast=paste(com,collapse="-"))
  write.table(res,file=gsub("pdf$","txt",strF),sep="\t",col.names = NA,quote=F)
  ## get the siginificant for each compared grp
  sigX <- !is.na(res$padj)&res$padj<padj&res$log2FoldChange< -logFC
  sigY <- !is.na(res$padj)&res$padj<padj&res$log2FoldChange> logFC
  ## plot-----------
  Col <- rep("gray",nrow(Data))
  Col[sigX] <- COL[2]
  Col[sigY] <- COL[1]
  x <- apply(Data[,pClass==com[2]],1,mean)
  y <- apply(Data[,pClass==com[1]],1,mean)
  
  pdf(strF,width=4,height=4,useDingbats = F)
  par(mar=c(2,2,2,0)+0.2,mgp=c(1,0,0),tcl=-0.02)
  ## plot scatter dots -----------
  imageCOL <- c("#FFFFFFFF","#BEBEBE","#B4B4B4","#AAAAAA","#A0A0A0","#969696",
                "#8C8C8C","#828282","#787878","#6E6E6E","#646464","#5A5A5A",
                "#505050","#464646","#3C3C3C","#323232","#282828","#1E1E1E",
                "#141414","#0A0A0A","#000000")
  plot(c(),c(),xlim=range(c(x,y)),ylim=range(c(x,y)),xlab=com[2],ylab=com[1],main=title)
  for(i in c("gray",COL)){
    if(i=="gray"){## plot not significant
      tryM <- try(f1 <- MASS::kde2d(x[Col==i],y[Col==i],n=50),silent=T)
      if(is.null(names(tryM))){
        index <- rep(T,sum(Col==i))
      }else{
        image(f1,col=imageCOL,add=T)
        imageZero <- diff(range(f1$z))/length(imageCOL)
        index <- apply(cbind(x[Col==i],y[Col==i]),1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
      }
    }else{
      index <- rep(T,sum(Col==i))
    }
    points(x[Col==i][index],y[Col==i][index],pch=20,col=i,cex=1)
  }
  mtext(paste(com[1],": ",sum(sigY),sep=""),3,line=-1,adj=0.1,col=COL[1])
  mtext(paste(com[2],": ",sum(sigX),sep=""),4,line=-1,adj=0.1,col=COL[2])
  ## plot sig genes -----------------
  if(sum(Col!="gray")>0){
    plot(c(),c(),xlim=range(c(x,y)),ylim=range(c(x,y)),xlab=com[2],ylab=com[1],main=title)
    for(i in COL)
      points(x[Col==i],y[Col==i],pch=20,col=paste(i,"80",sep=""),cex=1)
    text(x[Col!="gray"],y[Col!="gray"],rownames(Data)[Col!="gray"],cex=0.3)
  }
 
  ## volcano plots ----------------------
  s <- -log10(res$padj)
  s[is.infinite(s)] <- max(s[is.finite(s)])*1.05
  s[is.na(s)] <- min(s,na.rm=T)
  baseN <- 5*max(s)/pSize
  sLegend <- unique(10*round(max(s)/10)*c(0.25,0.5,0.75,1))
  logFC <- res$log2FoldChange
  meanT <- apply(cbind(x,y),1,mean)
  par(mar=c(2,2,1,7)+0.2,mgp=c(1,0,0),tcl=-0.02,xpd=T)
  plot(c(),c(),xlim=range(logFC),ylim=range(meanT),xlab=paste("logFC(",paste(com,collapse="/"),")"),ylab="Mean log2(TPM+1)",main=title)
  for(i in order(s,decreasing=T)) 
    if(s[i]>min(median(s),sLegend[1]))
      draw.circle(logFC[i],
                  meanT[i],
                  max(0.01,s[i]/baseN),
                  lty=0,
                  col=Col[i])
  ## legend
  loc <- par('usr')
  yLoc <- min(loc[3:4])+2*diff(loc[3:4])/3
  mtext(paste(com[1],": ",sum(sigY),sep=""),4,line=-1,adj=0.9,col=COL[1])
  mtext(paste(com[2],": ",sum(sigX),sep=""),2,line=-1,adj=0.9,col=COL[2])
  if(length(sLegend)>1){
    legendP <- legend("topright",legend=sLegend,yjust=1,title="-log10(p-value)",bty = "n",inset=c(-0.6,0),cex=0.8)#,pch=20,col="white"
    adjustX <- (legendP$rect$left-legendP$text$x)/2
    for(i in 1:length(sLegend)){
      draw.circle(legendP$text$x[i]+adjustX,legendP$text$y[i],sLegend[i]/baseN)
    }
  }
  ## mark the genes -----------------------
  if(sum(Col!="gray")>0){
    plot(c(),c(),xlim=range(logFC),ylim=range(meanT),xlab=paste("logFC(",paste(com,collapse="/"),")"),ylab="Mean log2(TPM+1)",main=title)
    for(i in order(s,decreasing=T)){
      if(Col[i]=="gray") next
      if(s[i]>min(median(s),sLegend[1])){
        draw.circle(logFC[i],
                    meanT[i],
                    max(0.01,s[i]/baseN),
                    lty=0,
                    col=paste(Col[i],"80",sep=""))
        text(logFC[i],meanT[i],rownames(Data)[i],cex=0.2)
      }
    }
  }
  ### -----
  a <- dev.off()
  return(rownames(res)[sigX|sigY])
}
DEG <- c()
for(i in unique(pClass)){
  for(j in unique(pClass)){
    if(i==j) next
    cat("\t\t",j,"vs",i,"\n")
    index <- pClass%in%c(i,j)
    DEG <- c(DEG,DESeq2.scatter(rawC[,index],pClass[index],opt$out,
                                COL[c(i,j)],cutOff=opt$minTPM,
                                showX=rawT[,index],
                                logFC=1,title="log2(TPM)"))
  }
}
## overall differential genes ------------------
DEG <- unique(DEG)
if(length(DEG)>3){
  cat("\n\tStep 3: summarize all differential genes to plot heatmap\n")
  heatCol <- c("#053061","#134C88","#2268AD","#3480B9","#4B98C5","#74B2D4",
               "#9BCAE0","#BDDAEA","#D8E8F1","#ECF2F5","#F8EFEA","#FBE0D1",
               "#FAC9B1","#F5AD8C","#E88B6E","#D96752","#C6413E","#B31B2C",
               "#8E0C25","#67001F")
  pdf(paste(opt$out,"/overall.heatmap.pdf",sep=""),height=9)
  pheatmap(log2(1+rawT[DEG,]),annotation_colors=list(grp=COL),fontsize_row=1,
           color = heatCol,
           annotation_col = data.frame(row.names=colnames(rawT),grp=pClass))
  pheatmap(log2(1+rawT[DEG,]),annotation_colors=list(grp=COL),scale="row",fontsize_row=1,
           color = heatCol,
           annotation_col = data.frame(row.names=colnames(rawT),grp=pClass))
  pheatmap(log2(1+rawT[DEG,]),annotation_colors=list(grp=COL),scale="row",cluster_cols=F,
           labels_row=rep("",length(DEG)),
           color = heatCol,
           annotation_col = data.frame(row.names=colnames(rawT),grp=pClass))
  a <- dev.off()
}

## ---------------
cat("\nDifferential analysis is done successfully!\n")







