#!/usr/bin/env Rscript
##################
## trimLen.R
##
##################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\n\tTrim fastqs for a given length.\n")
  cat("\tusage: trimLen.R /full/path/to/the/fastq/folder/ <length>\n")
  cat("\tlength: the number of leading bp to be kept\n")
  cat("\teg: trimLen.R /data/archive/z5ouyang/ 30")
  cat("\n\tResult: located in a new folder with the name appended with '_trim' plus the length after input folder name.\n")
  cat("\n\n")
  q()
}
strPath <- args[1]
strLength <- args[2]

strOut <- paste(dirname(strPath),"/",basename(strPath),"_trim",strLength,"/",sep="")
if(!dir.exists(strOut)) dir.create(strOut)
for(i in list.files(strPath,"fastq.gz$",full.names=T)){
  strF <- paste(strOut,gsub("fastq\\.gz","trim.fastq",basename(i)),sep="")
  system(paste("homerTools trim -len",strLength,i))
  system(paste("mv ",i,".trimmed ",strF,sep=""))
  system(paste("gzip",strF))
}
