#!/usr/bin/env Rscript
## initial ----
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse))
  stop("R packages of optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="atac", 
              help="The assay type atac, chip and (future: histone) [default= %default]", metavar="character"),
  make_option(c("-p", "--peak"), type="character", #default="", 
              help="The full path to the peak files separated by ','[default= %default]. Either '-p' or '-i' has to be provided", metavar="character"),
  make_option(c("-i", "--idr"), type="character", #default="", 
              help="The full path to a folder contains idr peaks of all groups defined sample file [default= %default]. Either '-p' or '-i' has to be provided", metavar="character"),
  make_option(c("-d", "--distal"), type="numeric", default="3000", 
              help="The distance from TSS to be considered enhancer region [default= %default]", metavar="character"),
  make_option(c("-t", "--track"), type="character", default="", 
              help="The name of the track to be generated. If not provided, no track will be generated. [default= %default]", metavar="character"),
  make_option(c("-c", "--homer"), type="character", default=" -size given -pc 3", 
              help="Additional commands for homer peak annotation (annotatePeaks.pl) [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="peakQuan.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/sample/definition/file is required.\n", call.=FALSE)
}
strSample <- args[1]
opt = parse_args(opt_parser,args[-1])
if(is.null(opt$peak) && is.null(opt$idr)){
  print_help(opt_parser)
  stop("Either '-p' or '-i' has to be provided.\n", call.=FALSE)
}
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
sInfo <- read.table(strSample,sep="\t",comment.char="",row.names=1,as.is = T)

## check the provided peaks----
strTmp <- paste(strOutput,"peakQuan_tmp/",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)
if(!is.null(opt$peak)){
  strPeaks <- unlist(strsplit(opt$peak,","))
}else{
  strPeaks <- paste(opt$idr,"/",rownames(sInfo),".idr",sep="")
}
cat("\n\tStep 1: check the provided peaks\n")
strPeak <- strPeaks
if(length(strPeaks)>1){
  strPeak <- paste(strTmp,"allPeaks",sep="")
  strCMD <- paste("mergePeaks",paste(strPeaks,collapse=" "),">",strPeak)
  cat("and two and more peak files are provided, merge them.\n")
  cat("\t\t",strCMD,"\n")
  system(paste(strCMD,"2>/dev/null"))
}
## annotate peaks on tag directories ------
cat("\n\tStep 2: Quantify the normalized and raw tag counts at the peaks\n")
strCMD <- paste("annotatePeaks.pl",strPeak,opt$genome,opt$homer,"-d",paste(unlist(strsplit(sInfo[,2],";")),collapse=" "),">",paste(strOutput,"allNormTags.txt",sep=""))
cat("\t\t",strCMD,"\n")
if(!file.exists(paste(strOutput,"allNormTags.txt",sep="")))
system(paste(strCMD,"2>/dev/null"))

strCMD <- paste("annotatePeaks.pl",strPeak,opt$genome,opt$homer,"-noadj -d",paste(unlist(strsplit(sInfo[,2],";")),collapse=" "),">",paste(strOutput,"allRawTags.txt",sep=""))
cat("\t\t",strCMD,"\n")
if(!file.exists(paste(strOutput,"allRawTags.txt",sep="")))
system(paste(strCMD,"2>/dev/null"))

## plot the peak distrubution in location for each sample -----------
cat("\n\tStep 3: plot the tag distrubution in peaks each sample based on allNormTags.txt\n")
sID <- c()
for(i in rownames(sInfo)) sID <- c(sID,paste(i,unlist(strsplit(sInfo[i,3],";")),sep="_"))
normTags <- read.table(paste(strOutput,"allNormTags.txt",sep=""),as.is=T,sep="\t",header=T,row.names=1,quote="",comment.char="")
normCounts <- normTags[,-(1:18)]
colnames(normCounts) <- sID
STATs <- matrix(0,nrow=3,ncol=length(sID),dimnames=list(c("Promoter",paste("Distal ",opt$distal/1000,"k",sep=""),"Not in Peaks"),sID))
STATs[1,] <- apply(normCounts[!is.na(normTags$Distance.to.TSS)&normTags$Distance.to.TSS<=opt$distal,],2,sum)/10000000
STATs[2,] <- apply(normCounts[is.na(normTags$Distance.to.TSS)|normTags$Distance.to.TSS>opt$distal,],2,sum)/10000000
STATs[3,] <- 1-apply(normCounts,2,sum)/10000000
COL <- c("#e41a1c","#377eb8","#4d4d4d")
pdf(paste(strOutput,"tagPeak_dist.pdf",sep=""),width=9)
par(mar=c(12,2,0,0)+0.2,mgp=c(0.5,0,0),tcl=-0.03)
plot(c(),c(),xlim=c(0,1),ylim=c(0,1),axes=F,xlab="",ylab="")
legend("center",rownames(STATs),fill=COL)
barplot(STATs,las=2,col=COL)
a <- dev.off()
## merge all tag directories of a group --------
cat("\n\tStep 4: merge tag directories for a group\n")
strMerge <- paste(strOutput,"mergeTag/",sep="")
if(!dir.exists(strMerge)) dir.create(strMerge)
conn <- file(paste(strMerge,"fragL.txt",sep=""),"w")
for(i in rownames(sInfo)){
  cat("\t\tmerge",i,"\n")
  strTag <- paste(strMerge,i,sep="")
  strCMD <- paste("makeTagDirectory",strTag,"-d",paste(unlist(strsplit(sInfo[i,2],";")),collapse=" "))
  cat("\t\t\t",strCMD,"\n")
  if(!dir.exists(strTag))
  {
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
    tagTmp <- read.table(paste(strTag,"/tagInfo_ori.txt",sep=""),sep="\t",as.is=T,header=T,check.names=F)
    tagTmp[grepl("fragmentLengthEstimate",tagTmp[,1]),1] <- paste("fragmentLengthEstimate",round(mean(fragL)),sep="=")
    write.table(tagTmp,file=paste(strTag,"/tagInfo.txt",sep=""),sep="\t",row.names=F,quote=F,na="")
  }
}
close(conn)
strTagDir <-paste(strMerge,rownames(sInfo),sep="")
## update the sequence depth to the total tags in peaks for ATAC -----
cat("\n\tStep 5: update the sequence depth to the total tags in peaks for ATAC only (escaped for chip)\n")
if(opt$assay=="atac"){
  strCMD <- paste("annotatePeaks.pl",strPeak,opt$genome,"-noann -noadj -d", paste(strTagDir,collapse=" "),">",paste(strMerge,"rawTags",sep=""))
  cat("\t\tannotate all merge tag directories on the peaks\n")
  cat("\t\t",strCMD,"\n")
  if(!file.exists(paste(strMerge,"rawTags",sep=""))){
    system(paste(strCMD,"2>/dev/null"))
    tags <- read.table(paste(strMerge,"rawTags",sep=""),header=T,sep="\t",quote="",comment.char="",check.names=F,row.names=1)
    colnames(tags) <- sapply(strsplit(colnames(tags)," "),head,1)
    write.table(t(t(apply(tags[,19:ncol(tags)],2,sum))),file=paste(strMerge,"totalNumTags",sep=""),sep="\t",col.names=F,quote=F)
    for(i in strTagDir){
      cat(i,"\n")
      tagInfo <- read.table(paste(i,"/tagInfo.txt",sep=""),sep="\t",as.is=T,header=T,check.names=F)
      tagInfo[1,3] <- sum(tags[,i])
      write.table(tagInfo,file=paste(i,"/tagInfo.txt",sep=""),sep="\t",row.names=F,quote=F,na="")
    }
  }

}
## make a hub of all samples -----------------
hubName <- opt$track#paste(Sys.info()['user'],gsub("\\.txt","",basename(strSample)),sep="_")
if(nchar(hubName)>2){
  cat("\n\tStep 6: create a hub",hubName,"for all merge tags\n")
  strCMD <- paste("makeMultiWigHub.pl",hubName,opt$genome,
                  "-color",paste(apply(col2rgb(sInfo[,1]),2,paste,collapse=","),collapse=" "),
                  "-force -d",paste(strTagDir,collapse=" "))
  cat("\t\t",strCMD,"\n")
  if(system(paste(strCMD,"2>/dev/null"))!=0) stop("Error in making a hub of all samples!")
  
  strHub <- paste("/homer_data/www/html/hubs/",hubName,"/",opt$genome,"/trackDb.txt",sep="")
  if(!file.exists(strHub)) stop("The Hub was NOT created successfully!\n")
  cat("\t\tAdd hub on ucsc genome browser:\n\n\n=================================================\n\t\thttp://homer.ucsd.edu/hubs/",
      hubName,"/hub.txt\n=================================================\n\n\n",sep="")
}
#file.copy(strHub,gsub("\\.txt","_ori.txt",strHub))
#hubs <- scan(strHub,what="character",sep="\n",quiet=T)
#for(i in rownames(sInfo)){
#  hubs[grep(paste("track",i),hubs)+6] <- paste("color",paste(as.vector(col2rgb(sInfo[i,1])), collapse = ","))
#}
#cat(paste(hubs,collapse="\n"),file=strHub)
## -----------------
cat("\nQuantification is done successfully!\n")




