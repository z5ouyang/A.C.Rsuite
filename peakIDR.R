#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cloud.r-project.org/")
if(!require(optparse))
  stop("R packages of reshape2, MASS, gplots or optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]),
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="atac",
              help="The assay type atac, chip and (future: histone) [default= %default]", metavar="character"),
  make_option(c("-c", "--homer"), type="character", default=" -L 0 -C 0 -fdr 0.9 -minDist 200 -size 200",
              help="Additional commands for homer peak calling [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/idr/sample/definition/file [options]",
                          option_list=option_list,prog="peakIDR.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/idr/sample/definition/file is required.\n", call.=FALSE)
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
strTmp <- paste(strOutput,"peakIDR_tmp",sep="")
if(!dir.exists(strTmp)) dir.create(strTmp)

strAllIDR <- c()
for(i in 1:nrow(exps)){
  grpID <- exps[i,1]
  cat("\nAnalysing ",grpID,"\n",sep="")
  tags <- unlist(strsplit(exps[i,3],";"))
  sIDs <- paste(grpID,unlist(strsplit(exps[i,4],";")),sep="_")
  message("\t\tTag directory number:",length(tags),"; sample name number:",length(sIDs))
  if(length(tags)!=length(sIDs)) stop("Tags are not matching with sample Name")

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
  }
  if(length(tags)<2){
    strAllIDR <- c(strAllIDR,paste(strOutput,grpID,".idr",sep=""))
    strCMD <- paste("cp",strPeak,tail(strAllIDR,1))
    cat("\n\tNO REPLICATES: using peak from",sIDs[1]," as IDR peaks\n")
    cat(strCMD)
    system(strCMD)
    next
  }
  ## pair-wised IDR
  strIDRs <- c()
  for(j in 1:(length(strPeak)-1)){
    for(k in (j+1):length(strPeak)){
      strIDRs <- c(strIDRs,paste(strTmp,"/",grpID,(1+length(strIDRs)),".idr",sep=""))
      strCMD <- paste("idr.R",paste(strPeak[c(j,k)],collapse = " "),"-p",tail(strIDRs,1))
      cat("\n\tStep 2: Finding pair-wised IDR peaks for",sIDs[j],"and",sIDs[k],"\n")
      cat("\t\t",strCMD,"\n")
      #if(!file.exists(tail(strIDRs,1)))
      system(paste(strCMD,"2>/dev/null"))
      if(!file.exists(tail(strIDRs,1))) stop(paste("IDR failed for",tail(strIDRs,1),"\nPlease check your peak correlation result or rerun the command again might solve the problem!"))
    }
  }
  ## merge all pair-wised IDR peaks
  strAllIDR <- c(strAllIDR,paste(strOutput,grpID,".idr",sep=""))
  strCMD <- paste("mergePeaks",paste(strIDRs,collapse = " "),">",tail(strAllIDR,1))
  cat("\n\tStep 3: Mergeing all pair-wised IDR peaks for",grpID,"\n")
  cat("\t\t",strCMD,"\n")
  #if(!file.exists(tail(strAllIDR,1)))
  system(paste(strCMD,"2>/dev/null"))
}

## merge all IDR peaks -----
cat("\n\tStep 4: Merging all IDR peaks\n")
strPeak <- paste(strOutput,gsub("\\.txt","",basename(strSample)),".idr",sep="")
strVenn <- paste(strOutput,gsub("\\.txt","",basename(strSample)),".venn",sep="")
strCMD <- paste("mergePeaks",paste(strAllIDR,collapse=" "),"-venn",strVenn,">",strPeak)
cat("\t\t",strCMD,"\n")
system(paste(strCMD,"2>/dev/null"))
cat("\n\tStep 5: plot venn for all IDR peaks\n")
strCMD <- paste("plotMergeVenn.R",strVenn)
cat("\t\t",strCMD,"\n")
system(paste(strCMD,"2>/dev/null"))

cat("\nThe IDR process are completed successfully!\n")
