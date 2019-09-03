#!/usr/bin/env Rscript
############################################
## fastqc.R
## 
#############################################

args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1){
  cat("\n\nO'Young:\n\tPerform the fastqc and generate report on fastq.gz files in a given folder.\n")
  cat("\tUsage: fastqc.R path/to/the/fastq.gz/folder\n")
  cat("\n\n")
  q()
}
if(nchar(system("which fastqc",intern=T))<6)
  stop("Please install fastqc first @ https://www.bioinformatics.babraham.ac.uk/projects/fastqc/")
if(nchar(system("which multiqc",intern=T))<6)
  stop("Please install multiqc first @ https://multiqc.info/docs/")

system(paste("fastqc ",args[1],"/*fastqc.gz",sep=""))
system(paste("multiqc ",args[1],"-o",args[1]))
cat(paste("\nQC on fastqs were run successfully!\nPlease double click ",args[1],"/multiqc_report.html for summary view!",sep=""))