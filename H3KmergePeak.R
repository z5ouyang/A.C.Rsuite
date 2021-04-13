#!/usr/bin/env Rscript
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cloud.r-project.org/")
if(!suppressWarnings(suppressMessages(require(ggplot2)))) install.packages("ggplot2",repos="https://cloud.r-project.org/")

args <- commandArgs(trailingOnly=TRUE)

option_list = list(
  make_option(c('--H3K','-k'), action="store", type="character", dest="H3K", default=NULL,
              help="The path to the differential peaks (from peakDiff.R) from broad peaks, such as H3K27ac at narrow peak center (such as ATAC) +/- 500", metavar="character"),
  make_option(c("--H3Kquan","-q"), action="store", type="character", dest="H3Kquan", default=NULL,
              help="The path to the narrow peak file which was used to quantify H3K broad signal, and also was used to merge with other narrow peaks", metavar="character"),
  make_option(c("--peak","-p"), action="store", type="character", dest="peak", default=NULL,
              help="The path to the merged narrow peaks (~200 from TF ChIP/ATAC), one of which is used to quantify the H3K27ac", metavar="character")
)
opt_parser = OptionParser("\n\t%prog [options]",
                          option_list=option_list,prog="H3KmergePeak.R")
if (length(args)<1){
  print_help(opt_parser)
  message("H3KmergePeak.R -h .../...act.peak -q .../...idr -p .../mergedPeaks...txt\n")
  q()
}
opt = parse_args(opt_parser,args)
if(!file.exists(opt$H3K)) stop(paste0("The peakDiff result file (",opt$H3K,") cannot be loacated!"))
if(!file.exists(opt$peak)) stop(paste0("The mergePeak file (",opt$peak,") cannot be loacated!"))

H3K <- read.table(opt$H3K,sep="\t",header=T,as.is=T,row.names=1)
mergeP <- read.table(opt$peak,sep="\t",as.is=T,row.names=1)
strF <- basename(opt$H3Kquan)
if(sum(grepl(strF,mergeP[,6]))==0) stop(paste0("The narrow peaks (",opt$H3Kquan,") used for H3K quantification is NOT used for mergePeak (",opt$peak,")!"))
mergeP <- mergeP[grepl(strF,mergeP[,6]),]
colIndex <- 7+ which(nchar(mergeP[mergeP[,7]==1,][1,8:ncol(mergeP)])>0)
if(sum(!rownames(H3K)%in%unlist(strsplit(mergeP[,colIndex],",")))!=0) stop(paste0("The narrow peaks (",opt$H3Kquan,") used for mergePeaks is NOT used for H3K quantification!"))
D <- data.frame(table(mergeP[sapply(mergeP[,colIndex],function(x)return(sum(unlist(strsplit(x,","))%in%rownames(H3K))>0)),6]))
colnames(D) <- c('Overlap','Freq')
write.csv(D,file=paste0(opt$H3K,".splitMerge.csv"),row.names=F,quote=F)
p <- ggplot(D, aes(x="", y=Freq, fill=Overlap)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+
  theme(legend.position='bottom',
        legend.title=element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"inches"),
        legend.text = element_text(size=5),
        legend.key.size=unit(0.2,'lines'))+
  guides(fill=guide_legend(ncol=1,override.aes = list(size=2)))
if(nrow(D)<10){
  p <- p+scale_fill_brewer(palette="Set1")
}else if(nrow(D)<13){
  p <- p+scale_fill_brewer(palette="Set3")
}
pdf(paste0(opt$H3K,".splitMerge.pdf"),width=3,height=4)
print(p)
a <- dev.off()
