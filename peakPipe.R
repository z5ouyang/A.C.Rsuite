#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse)) stop("R packages of optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-a", "--assay"), type="character", default="atac", 
              help="The assay type atac, chip and (future: histone) [default= %default]", metavar="character"),
  make_option(c("-c", "--homerCor"), type="character", default=" -minDist 200 -size 200", 
              help="Additional commands for homer peak calling for correlation[default= '%default']", metavar="character"),
  make_option(c("-i", "--homerIDR"), type="character", default=" -L 0 -C 0 -fdr 0.9 -minDist 200 -size 200", 
              help="Additional commands for homer peak calling [default= '%default']", metavar="character"),
  make_option(c("-q", "--homerQuan"), type="character", default=" -size given -pc 3", 
              help="Additional commands for homer peak annotation (annotatePeaks.pl) [default= '%default']", metavar="character"),
  make_option(c("-d", "--distal"), type="numeric", default="3000", 
              help="The distance from TSS to be considered enhancer region [default= %default]", metavar="numeric"),
  make_option(c("-f", "--cutoff"), type="numeric", default="4", 
              help="The cutoff value of tags of a peak [default= %default]", metavar="numeric"),
  make_option(c("-m", "--maxPeak"), type="numeric", default="10000", 
              help="The maximum number of peaks to include in the heatmap [default= %default]", metavar="character"),
  make_option(c("-l", "--logFC"), type="numeric", default="1", 
              help="The log fold change to be considered significant [default= %default]", metavar="numeric"),
  make_option(c("-t", "--track"), type="character", default=paste(Sys.info()['user'],gsub("\\.txt","",basename(args[1])),sep="_"), 
              help="The name of the track to be generated. If not provided, no track will be generated. [default= '%default']", metavar="character"),
  make_option(c("-p", "--padj"), type="numeric", default="0.05", 
              help="The significant p-adj to be considered significant [default= %default]", metavar="numeric")
)
##-----
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="peakPipe.R")
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

## obtain the alignment status --------------------
cat("\nGet the alignment status: --------------------\n")
strCMD <- paste("alignStats.R",strSample,">",paste(strOutput,"alignStats.txt",sep=""))
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in alignment status!")

## Calculate the pair-wised correlation --------------------
cat("\nCalculate the pair-wised correlation: --------------------\n")
strCMD <- paste("peakCor.R",strSample,"-o",paste(strOutput,"peakCor/",sep=""),"-g",opt$genome,"-a",opt$assay,paste("-c '",opt$homerCor,"'",sep=""))
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Calculating the pair-wised correlation!")
## Calculate the IDR peaks --------------------
cat("\nCalculate the IDR peaks for each group: --------------------\n")
strCMD <- paste("peakIDR.R",strSample,"-o",paste(strOutput,"/peakIDR/",sep=""),"-a",opt$assay,paste("-c '",opt$homerIDR,"'",sep=""))
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Calculating the IDR peaks for each group!")

## Quantify IDR peaks for each sample --------------------
cat("\nQuantify IDR peaks for each sample: --------------------\n")
strPeak <- paste(strOutput,"peakIDR/",gsub("\\.txt","",basename(strSample)),".idr",sep="")
strCMD <- paste("peakQuan.R",strSample,"-o",paste(strOutput,"peakQuan/",sep=""),"-g",opt$genome,"-a",opt$assay,
                "-p",strPeak,"-d",opt$distal,
                "-t",opt$track,
                paste("-c '",opt$homerQuan,"'",sep=""))
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Quantifying IDR peaks for each sample!")

## Find differential peaks -------------------
cat("\nFind differential peaks for pair-wised groups: --------------------\n")
strAnno <- paste(strOutput,"peakQuan/allRawTags.txt",sep="")
strCMD <- paste("peakDiff.R",strSample,
                "-o",paste(strOutput,"/peakDiff/",sep=""),"-g",opt$genome,"-a",opt$assay,
                "-q",strAnno,"-d",opt$distal,
                "-c",opt$cutoff,"-m",opt$maxPeak,
                "-l",opt$logFC,"-p",opt$padj)
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Finding differential peaks for pair-wised groups!")

warnings()
## finish -----------------
cat("\n",opt$assay,"peak analysis pipeline finished successfully!\n")















