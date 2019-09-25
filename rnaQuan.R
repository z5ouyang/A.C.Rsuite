#!/usr/bin/env Rscript
if(!require(optparse,quietly = T)) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!require(plot3D,quietly = T)) install.packages("plot3D",repos="https://cran.cnr.berkeley.edu/")
if(!require(pheatmap,quietly = T)) install.packages("pheatmap",repos="https://cran.cnr.berkeley.edu/")
if(!require(MASS,quietly = T)) install.packages("MASS",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse) || !require(plot3D) || !require(pheatmap) || !require(MASS,quietly = T))
  stop("R packages of optparse, pheatmap, MASS or plot3D cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-l", "--length"), type="numeric", default="200", 
              help="The gene length cut-off [default= %default]", metavar="character"),
  make_option(c("-t", "--track"), type="character", default="", 
              help="The name of the track to be generated. If not provided, no track will be generated. [default= %default]", metavar="character"),
  make_option(c("-c", "--homer"), type="character", default="", 
              help="Additional commands besides '-condenseGenes -count exons -noadj' for homer annotation calling [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="rnaQuan.R")
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
## processing ------
## table 4 columns
## Column 1: grpname
## Column 2: col(rgb)
## Column 3: sampleTag1;sampleTag2;...
## Column 4: sID1;sID2;...
## obtain raw counts -----------
cat("\n\tStep 1: obtain the condenseGenes raw counts\n")
sInfo <- read.table(strSample,comment.char="",as.is=T,sep="\t",row.names=1)#
#print(sInfo)
strCMD <- paste("analyzeRepeats.pl rna",opt$genome,"-condenseGenes -count exons -noadj",opt$homer,"-d",paste(gsub(";"," ",sInfo[,2]),collapse = " "),">",paste(strOutput,"HOMER.rawCount.txt",sep=""))
cat("\t\t",strCMD,"\n")
system(paste(strCMD,"2>/dev/null"))

## obtain TPM counts -----------
cat("\n\tStep 2: obtain the TPM without condenseGenes then selecting transcripts according to raw counts\n")
strCMD <- paste("analyzeRepeats.pl rna",opt$genome,"-count exons -tpm",opt$homer,"-d",paste(gsub(";"," ",sInfo[,2]),collapse = " "),">",paste(strOutput,"HOMER.rawTPM.txt",sep=""))
cat("\t\t",strCMD,"\n")
system(paste(strCMD,"2>/dev/null"))
## matching TPM with condenseGene from raw counts
rawC <- read.table(paste(strOutput,"HOMER.rawCount.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
rawT <- read.table(paste(strOutput,"HOMER.rawTPM.txt",sep=""),sep="\t",as.is=T,row.names=1,check.names = F,comment.char = "",quote="",header=T)
rawT <- rawT[rownames(rawC),]
write.table(rawT,file=paste(strOutput,"HOMER.rawTPM.txt",sep=""),sep="\t",col.names = NA,quote=F)

## remove the short genes ---------------
cat("\n\tStep 3: remove the short genes and annotation information\n")
COL <- sInfo[,1]
names(COL) <- rownames(sInfo)
names(sID) <- sID <- colnames(rawC)[1:7]
pClass <- c()
for(i in rownames(sInfo)){
  one <- paste(i,unlist(strsplit(sInfo[i,3],";")),sep="_")
  names(one) <- unlist(strsplit(sInfo[i,2],";"))
  pClass <- c(pClass,rep(i,length(one)))
  sID <- c(sID,one)
}
dimnames(rawC) <- list(sapply(strsplit(rawC[,7],"\\|"),head,1),
                       sID[sapply(strsplit(colnames(rawC)," "),head,1)])
dimnames(rawT) <- list(sapply(strsplit(rawT[,7],"\\|"),head,1),
                       sID[sapply(strsplit(colnames(rawT)," "),head,1)])

if(is.data.frame(rawC)) rawC <- as.matrix(rawC[,c(5,8:ncol(rawC))])
if(is.data.frame(rawT)) rawT <- as.matrix(rawT[,c(5,8:ncol(rawT))])
## Plotting the gene length vs expression
cutLength <- opt$length
strMain <- gsub("\\.txt","",basename(strSample))
cat("\t\tplot the gene length vs expression, and length cutoff=",cutLength,"\n")
plotGeneLength <- function(X,funs,strM,cutLength){
  for(i in names(funs)){
    y <- apply(X[,-1],1,funs[[i]])
    plot(X[y>0,"Length"],y[y>0],xlab="gene length",ylab=paste(strM,i),type="p",pch=20,log="xy")
    lines(rep(cutLength,2),c(1e-2,1e5),col=2)
    axis(1,cutLength,paste(cutLength),col=2,col.axis=2)
  }
}
pdf(paste(strOutput,strMain,".geneLength.pdf",sep=""),width=9)
par(mar=c(5,5,0.5,0.5)+0.1)
if(sum(colnames(rawC)=="Length")>0){
  plotGeneLength(rawC,list(mean=mean,median=median,SD=sd),"raw counts",cutLength)
}
if(sum(colnames(rawT)=="Length")>0){
  plotGeneLength(rawT,list(mean=mean,median=median,SD=sd),"TPM",cutLength)
}
a <- dev.off()
rawC <- rawC[rawC[,1]>cutLength,colnames(rawC)!="Length"]
rawT <- rawT[rawT[,1]>cutLength,colnames(rawT)!="Length"]
sID <- colnames(rawC)
write.table(rawT,file=paste(strOutput,"rawT.txt",sep=""),sep="\t",col.names = NA,quote=F)
write.table(rawC,file=paste(strOutput,"rawC.txt",sep=""),sep="\t",col.names = NA,quote=F)

## plot PCA of TPM on all genes of all samples -----------------
cat("\n\tStep 4: plot PCA of TPM on all genes of all samples\n")
pCol <- COL[pClass]
pPch <- rep(20,length(pCol))
main <- paste(gsub("\\txt","",basename(strSample)),nrow(rawT),"genes")
logT <- log2(1+rawT)
logT <- logT[apply(logT,1,sd)>0,]
Data.pca <- summary(prcomp(t(logT),scale=TRUE,retx=TRUE))
pca.x <- Data.pca$x[,1]
pca.y <- Data.pca$x[,2]
pca.z <- Data.pca$x[,3]
pdf(paste(strOutput,"overall.PCA.pdf",sep=""))
plot(pca.x,pca.y,
     xlab=paste("PC1 (v:",round(Data.pca$importance[2,1]*100,1),"%)",sep=""),
     ylab=paste("PC2 (v:",round(Data.pca$importance[2,2]*100,1),"%)",sep=""),
     main=main,col=pCol,pch=pPch)#pColor
plot(pca.x,pca.y,
     xlab=paste("PC1 (v:",round(Data.pca$importance[2,1]*100,1),"%)",sep=""),
     ylab=paste("PC2 (v:",round(Data.pca$importance[2,2]*100,1),"%)",sep=""),
     main=main,col=pCol,pch=pPch)#pColor
text(pca.x,pca.y+diff(range(pca.y))/40,rownames(Data.pca$x),col=pCol,cex=0.3)

theta <- seq(0,360,length.out=36)
phi <- rep(45,length(theta))
par(col.axis="gray")
for(i in 1:length(theta)){
  #perspbox(pca.x,pca.y,pca.z,col.axis="gray",phi=phi[i],theta=theta[i])
  scatter3D(pca.x,pca.y,pca.z,bty="u",col.axis="gray90",col.grid="gray95",cex=2,
            colvar=NULL,#add=T,
            col=pCol,pch=pPch,
            main=main,phi=phi[i],theta=theta[i],
            xlab=paste("PC1 (v:",round(Data.pca$importance[2,1]*100,1),"%)",sep=""),
            ylab=paste("PC2 (v:",round(Data.pca$importance[2,2]*100,1),"%)",sep=""),
            zlab=paste("PC3 (v:",round(Data.pca$importance[2,3]*100,1),"%)",sep="")
  )
  text3D(pca.x,pca.y,pca.z,
         rownames(Data.pca$x),
         col=pCol,cex=0.3,add=T)
  
}
a <- dev.off()
## Calculate the correlations ---------------------------
cat("\n\tStep 5: calculate the correlations of TPM among replicates\n")
heatCol <- c("#99000D","#A80712","#B80F17","#C8161C","#D42120","#DF2C25",
             "#EB372A","#F14432","#F5533B","#F96245","#FB7050","#FB7C5C",
             "#FB8969","#FC9676","#FCA385","#FCB094","#FCBDA3","#FCCAB5",
             "#FDD7C7","#FEE5D9")
corT <- cor(logT,logT)
pdf(paste(strOutput,strMain,".sampleCor.pdf",sep=""),width=9,height=9)
pheatmap(corT,rev(heatCol),cellwidth=20,cellheight=20,
         cluster_rows=F,cluster_cols=F,
         fontsize=8,
         annotation_row=data.frame(pClass=pClass,row.names=rownames(corT)),
         annotation_col=data.frame(pClass=pClass,row.names=colnames(corT)),
         annotation_colors=list(pClass=COL),
         display_numbers=T,number_color="black")
## replicates for each group
cutOff <- 0
limFun <- function(x){
  step <- diff(range(x))/20
  return(range(x)-c(step,-step))
}
imageCOL <- c("#FFFFFFFF","#3300FF","#2D1CFF","#2839FF","#2255FF","#1C71FF","#178EFF","#11AAFF",
              "#0BC6FF","#06E3FF","#00FFFF","#00FFFF","#17FFE3","#2DFFC6","#44FFAA","#5BFF8E",
              "#71FF71","#88FF55","#9FFF39","#B5FF1C","#CCFF00",
              "#CCFF00","#D2F400","#D7E800","#DDDD00","#E3D200","#E8C600","#EEBB00","#F4B000",
              "#F9A400","#FF9900","#FF9900","#F68800","#EC7700","#E36600","#D95500","#D04400",
              "#C63300","#BD2200","#B31100","#AA0000")
X <- logT
for(type in unique(pClass)){
  cat("\t\tworking on",type,"\n")
  ids <- colnames(logT)[pClass==type]
  COR <- matrix(1,nrow=length(ids),ncol=length(ids),dimnames=list(ids,ids))
  for(i in ids){
    for(j in ids){
      if(grep(paste(i,"$",sep=""),ids)>=grep(paste(j,"$",sep=""),ids)) next
      index <- apply(X[,c(i,j)],1,max)>cutOff
      tmp <- X[index,c(i,j)]
      COR[i,j] <- COR[j,i] <- cor(tmp[,1],tmp[,2])
      xlim <- limFun(tmp[,1])
      ylim <- limFun(tmp[,2])
      tryM <- try(f1 <- kde2d(tmp[,1],tmp[,2],n=200),silent=T)
      if(is.null(names(tryM))){
        plot(c(),c(),xlab=i,ylab=j,xlim=xlim,ylim=ylim)
        index <- rep(T,nrow(tmp))
      }else{
        image(f1,col=imageCOL,xlab=i,ylab=j,xlim=xlim,ylim=ylim)
        imageZero <- diff(range(f1$z))/length(imageCOL)
        index <- apply(tmp,1,function(x,fit,cutZero){return(fit$z[sum((x[1]-fit$x)>=0),sum((x[2]-fit$y)>=0)]<cutZero)},f1,imageZero)
      }
      points(tmp[index,1],tmp[index,2],col=imageCOL[2],pch=20)
      lines(range(c(xlim,ylim)),range(c(xlim,ylim)),col="gray")
      legend("topleft",paste("r=",round(COR[i,j],3)),box.lty=0)
    }
  }
  pheatmap(COR,rev(heatCol),cellwidth=20,cellheight=20,
           cluster_rows=F,cluster_cols=F,
           fontsize=8,
           display_numbers=T,number_color="black")
  
}
a <- dev.off()

## Merge all group tag directories -------
cat("\tStep 6: merge all tag directories for each group\n")
strMerge <- paste(strOutput,"mergeTag/",sep="")
if(!dir.exists(strMerge)) dir.create(strMerge)
conn <- file(paste(strMerge,"fragL.txt",sep=""),"w")
for(i in rownames(sInfo)){
  cat("\t\tmerge",i,"\n")
  strTag <- paste(strMerge,i,sep="")
  strCMD <- paste("makeTagDirectory",strTag,"-d",paste(unlist(strsplit(sInfo[i,2],";")),collapse=" "))
  cat("\t\t\t",strCMD,"\n")
  system(paste(strCMD,"2>/dev/null"))
  ## get average fragment length
  cat("\t\t\tUpdate the average fragment length to the merged tag directory\n")
  fragL <- c()
  for(j in unlist(strsplit(sInfo[i,2],";"))){
    res <- read.table(paste(j,"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T)
    fragL <- c(fragL,as.numeric(gsub("fragmentLengthEstimate=","",res[grepl("fragmentLengthEstimate",res[,1]),1])))
    cat(basename(j),"\t",tail(fragL,1),"\n",file=conn)
  }
  system(paste("cp ",strTag,"/tagInfo.txt ",strTag,"/tagInfo_ori.txt",sep=""))
  tagInfo <- read.table(paste(strTag,"/tagInfo_ori.txt",sep=""),sep="\t",as.is=T,header=T,check.names=F)
  tagInfo[grepl("fragmentLengthEstimate",tagInfo[,1]),1] <- paste("fragmentLengthEstimate",round(mean(fragL)),sep="=")
  write.table(tagInfo,file=paste(strTag,"/tagInfo.txt",sep=""),sep="\t",row.names=F,quote=F,na="")
}
close(conn)
strTagDir <-paste(strMerge,rownames(sInfo),sep="")
## make a hub of all samples -----------------
hubName <- opt$track
if(nchar(hubName)>2){
  cat("\tStep 7: create a hub",hubName,"for all merge tags, will try to overwrite if it exists\n")
  strCMD <- paste("makeMultiWigHub.pl",hubName,opt$genome,
                  "-force -d",paste(strTagDir,collapse=" "))
  cat("\t\t",strCMD,"\n")
  system(paste(strCMD,"2>/dev/null"))
  strHub <- paste("/homer_data/www/html/hubs/",hubName,"/",opt$genome,"/trackDb.txt",sep="")
  if(!file.exists(strHub)) stop("The Hub was NOT created successfully!\n")
  hubs <- scan(strHub,what="character",sep="\n",quiet=T)
  for(i in rownames(sInfo)){
    hubs[grep(paste("track",i),A)+6] <- paste("color",paste(as.vector(col2rgb(sInfo[i,1])), collapse = ","))
  }
  cat(paste(hubs,collapse="\n"),file=strHub)
  cat("\t\tAdd hub on ucsc genome browser:\n\t\thttp://homer.ucsd.edu/hubs/",hubName,"/hub.txt\n",sep="")
}

## -------------------
cat("\nQuantification is done successfully!\n")






