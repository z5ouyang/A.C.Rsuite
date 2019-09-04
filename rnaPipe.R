#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!require(optparse)) stop("R packages of optparse cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)
option_list = list(
  make_option(c("-o", "--out"), type="character", default=dirname(args[1]), 
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="mm10", 
              help="The genome of the sequeces came from, mm10, hg38 [default= %default]", metavar="character"),
  make_option(c("-l", "--length"), type="numeric", default="200", 
              help="The gene length cut-off [default= %default]", metavar="character"),
  make_option(c("-t", "--track"), type="character", default=paste(Sys.info()['user'],gsub("\\.txt","",basename(args[1])),sep="_"), 
              help="The name of the track to be generated. If empty string ('') provided, no track will be generated. [default= '%default']", metavar="character"),
  make_option(c("-m", "--minTPM"), type="numeric", default="8", 
              help="The minumal TPM required in at least 2 samples", metavar="character"),
  make_option(c("-c", "--homer"), type="character", #default=" ", 
              help="Additional commands besides '-condenseGenes -count exons -noadj' for homer annotation calling [default= %default]", metavar="character")
)
opt_parser = OptionParser("\n\t%prog path/to/the/sample/definition/file [options]",
                          option_list=option_list,prog="rnaPipe.R")
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
## obtain the alignment status --------------------
cat("\nGet the alignment status: --------------------\n")
strCMD <- paste("alignStats.R",strSample,">",paste(strOutput,"alignStats.txt",sep=""))
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in alignment status!")

## quantify the genes --------------------
cat("\nQuantification of genes: --------------------\n")
strCMD <- paste("rnaQuan.R",strSample,"-o",paste(strOutput,"rnaQuan/",sep=""),"-g",opt$genome,"-l",opt$length,"-t",opt$track)
if(!is.null(opt$homer))
  strCMD <- paste("rnaQuan.R",strSample,"-o",paste(strOutput,"rnaQuan/",sep=""),"-g",opt$genome,"-l",opt$length,"-t",opt$track,"-c",opt$homer)
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Quantifying genes for each sample!")

## Differential analysis ------------------
cat("\nDifferential analysis: --------------------\n")
strCMD <- paste("rnaDiff.R",strSample,"-o",paste(strOutput,"rnaDiff/",sep=""),
                "-c",paste(strOutput,"rnaQuan/rawC.txt",sep=""),
                "-t",paste(strOutput,"rnaQuan/rawT.txt",sep=""),
                "-m",opt$minTPM)
cat(strCMD,"\n")
if(system(strCMD)!=0) stop("Error in Finding differential expressed genes for pair-wised groups!")

## finish -----------------
cat("\nRNAseq pipeline finished successfully!\n")




