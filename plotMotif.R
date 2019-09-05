#!/usr/bin/env Rscript
#########################################
## plotMotif.R
##
#########################################
## findMotifsGenome.pl XXX_induced.txt genome ./XXX/ -size given -mknown /home/z5ouyang/src/data/all_threshold_0.5.motif -bg ./XXX_bg.txt
##
if(!suppressWarnings(suppressMessages(require(optparse)))) install.packages("optparse",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(grImport,quietly=T)))) install.packages("grImport",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(gridExtra,quietly=T)))) install.packages("gridExtra",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(pheatmap,quietly=T)))) install.packages("pheatmap",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(RColorBrewer,quietly=T)))) install.packages("RColorBrewer",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(colorspace,quietly=T)))) install.packages("colorspace",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(ggplot2,quietly=T)))) install.packages("ggplot2",repos="https://cran.cnr.berkeley.edu/")
if(!suppressWarnings(suppressMessages(require(htmltab,quietly = T)))) install.packages("htmltab",repos="https://cran.cnr.berkeley.edu/")

if(!require(optparse)||!require(grImport)||!require(gridExtra)||!require(pheatmap)||!require(RColorBrewer)||!require(colorspace)||!require(ggplot2)||!require(htmltab))
  stop("R packages of optparse, grImport, gridExtra, pheatmap, RColorBrewer, colorspace, ggplot2 or htmltab cannot be installed!")

args <- commandArgs(trailingOnly=TRUE)
option_list = list(
)
opt_parser = OptionParser("\n\t%prog path/to/the/motif/file\nMotif file is a tab separated file with 4 columns and header:\n
                          \theader:motif\tcol\tknown\tdenovo
                          \t1: path to the motif folder;\n
                          \t2: color of this motif set;\n
                          \t3: index of known motif separated by ',' at least 3\n
                          \t4: index of de novo motif separated by ',' at least 3",
                          option_list=option_list,prog="peakQuan.R")
if (length(args)<1){
  print_help(opt_parser)
  stop("path/to/the/motif/file\nMotif file is required.\n", call.=FALSE)
}

## input ----
motifchar <- 200
strInput <- args[1]
res <- read.table(strInput,as.is=T,sep="\t",comment.char = "",header=T)
strDir <- res[,1]
COL <- res[,2]
kIndex <- strsplit(res[,3],",")
hIndex <- strsplit(res[,4],",")
if(length(kIndex)<3||length(hIndex)<3) stop("Please provide at least 3  indexes for motif!")
names(kIndex) <- names(hIndex) <- basename(strDir)
strPDF <- paste(strInput,"_allMotif.pdf",sep="")

pdf(strPDF,width=9)
par(bg=NA)
## known motif -------
selMotif <- motifP <- motifR <- c()
logos <- list()
#cat(strDir)
for(i in strDir){
  i <- gsub("~",normalizePath("~"),i)
  strMotif <- paste(i,"/knownResults.txt",sep="")
  if(file.exists(strMotif)){
    cat("Plotting known motif table for",basename(i),"\n")
    one <- read.table(strMotif,sep="\t",header=T,as.is=T,check.names=F,comment.char="")
    one <- one[!duplicated(one[,1]),c(1,4,7,9),drop=F]
    dimnames(one) <- list(one[,1],c("motif",paste(basename(i),c("logP","target","bg"),sep="_")))
    ## selected known motifs and plot them as a table
    selMotif <- c(selMotif,one[as.numeric(kIndex[[basename(i)]]),1])
    oneTable <- list(textGrob("Name"),textGrob("Motif"),textGrob("-logP"),textGrob("Target(%)/BG(%)"))
    lay <- matrix(c(1,1,2,2,2,3,4),nrow=1)
    for(j in as.numeric(kIndex[[basename(i)]])){
      strLogo <- paste(i,"/knownResults/known",j,".logo.svg",sep="")
      if(!file.exists(gsub("svg$","ps",strLogo))){
        #cat("generating",strLogo,"\n")
        #cat(paste("inkscape ",strLogo," --export-ps=",gsub("svg$","ps",strLogo),sep=""),"\n")
        system(paste("inkscape ",strLogo," --export-ps=",gsub("svg$","ps",strLogo),sep=""))
      }
      if(!file.exists(gsub("svg$","xml",strLogo))) PostScriptTrace(gsub("svg$","ps",strLogo),gsub("svg$","xml",strLogo))
      logos[[sapply(strsplit(one[j,1],"\\/"),head,1)]] <- readPicture(gsub("svg$","xml",strLogo))
      oneTable <- c(oneTable,
                    list(textGrob(sapply(strsplit(one[j,1],"\\/"),head,1))),
                    list(pictureGrob(logos[[sapply(strsplit(one[j,1],"\\/"),head,1)]])),
                    list(textGrob(-1*one[j,2])),
                    list(textGrob(gsub("%","",paste(one[j,3],one[j,4],sep="/")))))
      lay <- rbind(lay,max(lay)+lay[1,])
    }
    if(nrow(lay)<15) lay <- rbind(lay,matrix(1+max(lay),nrow=15-nrow(lay),ncol=ncol(lay)))
    #save(logos,one,oneTable,kIndex,file="t.RData")
    plot.new()
    grid.draw(arrangeGrob(grobs=oneTable,layout_matrix=lay,top=basename(i)))#ncol=4
    #stop()
    # logP
    motifP <- merge(motifP,-one[,2,drop=F],by="row.names",all=T)
    rownames(motifP) <- motifP[,1]
    motifP <- motifP[,-1,drop=F]
    # enrichment
    one <- one[,3:4]
    motifR <- merge(motifR,matrix(as.numeric(gsub("%","",as.matrix(one))),ncol=ncol(one),dimnames=dimnames(one)),
                    by="row.names",all=T,sort=F)
    rownames(motifR) <- motifR[,1]
    motifR <- motifR[,-1,drop=F]
  }
}
motifP <- motifP[unique(selMotif),,drop=F]
motifR <- motifR[unique(selMotif),,drop=F]

# plot logP heatmap
cat("Plotting known motif heatmap\n")
if(length(motifP)==0) stop("Cannot locate any motif analyses result")
motifP[is.na(motifP)] <- 0
rownames(motifP) <- substr(rownames(motifP),1,apply(cbind(nchar(rownames(motifP)),motifchar),1,min))
if(is.null(COL)){
  if(ncol(motifP)<=8){
    COL <- brewer.pal(n = 8, name ="Dark2")[1:ncol(motifP)]
  }else{
    COL <-rainbow_hcl(ncol(motifP))
  }
}
#COL <- brewer.pal(n = 8, name ="Dark2")[1:ncol(motifP)]
names(COL) <- colnames(motifP)
#print(COL)
clusterCall <- function(hc, mat){
  sv = mat[,ncol(mat)]-mat[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
sMotif <- sort(motifP[motifP!=0])
heatCOL <- colorRampPalette(brewer.pal(n = 7, name ="Oranges"))(min(length(sMotif)-2,10))#RdYlBu
br <- c(0,sMotif[floor(seq(1,length(sMotif)-1,length.out=length(heatCOL)-1))]-min(diff(sMotif))/10,max(sMotif)+min(diff(sMotif))/10)
pHeat <- pheatmap(motifP,cluster_cols=F,annotation_colors=list(grp=COL),clustering_method="ward.D",
                  color = heatCOL,clustering_callback = clusterCall,
                  breaks=br,
                  annotation_col=data.frame(row.names=colnames(motifP),
                                            grp=colnames(motifP)))#,silent = T
#save(pHeat,logos,file="pHeatLogo.RData")
# plot known motif logo
gT <- pHeat$gtable
tmp <- list()
lay <- c()
#cat(paste(sapply(strsplit(gT$grobs[[4]]$label,"\\/"),head,1),collapse="\n"))
for(i in sapply(strsplit(gT$grobs[[4]]$label,"\\/"),head,1)){
  #cat("logo",i,"\n")
  tmp <- c(tmp,
           list(textGrob(i,x = unit(0, "npc"),just="left")),
           list(pictureGrob(logos[[i]])))
  lay <- rbind(lay,max(c(0,lay))+c(1,2,2,2))
}
#save(gT,file="gT.RData")
gT$grobs[[4]] <- arrangeGrob(grobs=tmp,layout_matrix=lay)#,ncol=2
plot.new()
grid.draw(gT)
#save(motifR,motifP,file="t.RData")
# plot bubble plot for enrichment
cat("Plotting known motif bubble plots\n")
X <- data.frame()
mNames <- sapply(strsplit(rownames(motifP),"\\/"),head,1)
if(sum(duplicated(mNames))>0){
  cat("duplicated motif names:\n",paste(rownames(motifP)[mNames%in%mNames[duplicated(mNames)]],collapse="\n"),"\n\n")
  mNames[duplicated(mNames)] <- paste(mNames[duplicated(mNames)],"_1",sep="")
  cat("Skip bubble plot ...\n")
}else{
  rownames(motifP) <- mNames
  for(i in gsub("_bg$","",grep("_bg$",colnames(motifR),value=T))){
    #print(paste(i,"bg",sep="_"))
    X <- rbind(X,data.frame(grp=i,#motif=substr(rownames(motifR),1,apply(cbind(nchar(rownames(motifR)),motifchar),1,min)),
                            motif = factor(sapply(strsplit(rownames(motifR),"\\/"),head,1),
                                           levels=rev(sapply(strsplit(pHeat$gtable$grobs[[4]]$label,"\\/"),head,1))),
                            enrichment=motifR[,paste(i,"target",sep="_")]/motifR[,paste(i,"bg",sep="_")],
                            target=motifR[,paste(i,"target",sep="_")],
                            sig=motifP[sapply(strsplit(rownames(motifR),"\\/"),head,1),paste(i,"logP",sep="_")]))
  }
  
  print(ggplot(X,aes(x=grp,y=motif))+geom_point(aes(size=enrichment,colour=target))+scale_size_continuous(range = c(1,10))+
          scale_color_gradient(low="#fee5d9", high="#a50f15")+
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line =element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1,size=15),
                axis.text.y = element_text(hjust = 1,size=15),
                panel.background = element_blank()))
  
  print(ggplot(X,aes(x=grp,y=motif))+geom_point(aes(size=enrichment,colour=sig))+scale_size_continuous(range = c(1,10))+
          scale_color_gradient(low=heatCOL[1], high=tail(heatCOL,1))+
          theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line =element_blank(),
                axis.text.x = element_text(angle = 45, hjust = 1,size=15),
                axis.text.y = element_text(hjust = 1,size=15),
                panel.background = element_blank()))
}

## homer motif -------
#ori <- par(mar=c(2,7,max(c(0,(30-length(deNovo)*2))),1)+0.1,mgp=c(0.5,0,0),tcl=-0.1)
names(COL) <- gsub("_logP","",names(COL))
#print(COL)
#print(strDir)
for(i in strDir){
  i <- gsub("~",normalizePath("~"),i)
  deNovo <- as.numeric(hIndex[[basename(i)]])
  if(length(deNovo)<=0) next
  ori <- par(mar=c(2,7,max(c(0,(30-length(deNovo)*2))),1)+0.1,mgp=c(0.5,0,0),tcl=-0.1)
  strMotif <- paste(i,"/homerResults.html",sep="")
  if(!file.exists(strMotif)) next
  cat("Plotting deNovo motifs for",basename(i),"\n")
  motifs <- htmltab(strMotif,1,1)
  # obtain the motif names
  motifID <- c()
  oneTable <- list(textGrob("Motif"),textGrob("-logP"),textGrob("Target(%)/BG(%)"),textGrob("Matches"))
  lay <- matrix(c(1,1,1,2,3,4,4,4),nrow=1)
  for(j in deNovo){
    res <- scan(paste(i,"/homerResults/motif",j,".info.html",sep=""),character(),quiet=T)
    motifID <- c(motifID,
                 paste(gsub("<H4>","",sapply(strsplit(head(res[grep("<H4>",res)],3),"\\/"),head,1)),collapse="/"))
    
    strLogo <- paste(i,"/homerResults/motif",j,".logo.svg",sep="")
    system(paste("inkscape ",strLogo," --export-ps=",gsub("svg$","ps",strLogo),sep=""))
    PostScriptTrace(gsub("svg$","ps",strLogo),gsub("svg$","xml",strLogo))
    logo <- readPicture(gsub("svg$","xml",strLogo))
    oneTable <- c(oneTable,
                  list(pictureGrob(logo)),
                  list(textGrob(-1*as.numeric(motifs[j,"log P-pvalue"]))),
                  list(textGrob(gsub("%","",paste(motifs[j,"% of Targets"],motifs[j,"% of Background"],sep="/")))),
                  list(textGrob(tail(motifID,1),x=unit(0, "npc"),just="left")))
    lay <- rbind(lay,max(lay)+lay[1,])
  }
  if(nrow(lay)<15) lay <- rbind(lay,matrix(1+max(lay),nrow=15-nrow(lay),ncol=ncol(lay)))
  plot.new()
  grid.draw(arrangeGrob(grobs=oneTable,layout_matrix=lay,top=basename(i)))#ncol=4

  # obtain the motif info
  motifs <- motifs[rev(deNovo),]
  logP <- -as.numeric(motifs[,"log P-pvalue"])
  maxP <- ceiling(max(logP)/10)*10
  maxR <- maxP/3
  stepN <- 4
  maxTarget <- stepN*ceiling(max(as.numeric(gsub("%","",motifs[,"% of Targets"])))/stepN)
  # significance p-value bar plot
  y <- barplot(logP,horiz=T,xlab="-log(p-value)",las=1,main=basename(i),xlim=c(0,maxP+maxR),col=COL[basename(i)])
  #text(rep(maxP/2,length(y)),y,rev(motifID),col="gray50")
  a <- 255-mean(as.vector(col2rgb(COL[basename(i)])))
  text(rep(0,length(y)),y,rev(motifID),pos=4,col=rgb(a,a,a,maxColorValue=255))
  # target/bg ratio line plot
  Col <- c(frame="black",tg="#006d2c",bg="#74c476")
  axis(3,maxP+maxR*c(0:stepN)/stepN,paste(c(0:stepN)*maxTarget/stepN,"%",sep=""),col=Col["frame"],col.axis=Col["frame"])
  #axis(3,c(maxP,maxP+maxR/4,maxP+maxR/2,maxP+maxR*3/4,maxP+maxR),c("0%","25%","50%","75%","100%"),col=Col["frame"],col.axis=Col["frame"])
  mtext("Percentage in sequence",3,1,at=maxP+maxR,col=Col["frame"],adj=1)
  lines(as.numeric(gsub("%","",motifs[,"% of Targets"]))*maxR/maxTarget+maxP,y,col=Col["tg"])
  points(as.numeric(gsub("%","",motifs[,"% of Targets"]))*maxR/maxTarget+maxP,y,col=Col["tg"],pch=15)
  lines(as.numeric(gsub("%","",motifs[,"% of Background"]))*maxR/maxTarget+maxP,y,col=Col["bg"],lty=2)
  points(as.numeric(gsub("%","",motifs[,"% of Background"]))*maxR/maxTarget+maxP,y,col=Col["bg"],pch=16)
  lines(rep(maxP,2),c(0,max(y)*2),col=Col["frame"],lty=2)
  legend("bottomright",names(Col)[-1],col=Col[-1],lty=1:2,pch=15:16,box.lty=0,text.col=Col[-1])
}
par(ori)
a <- dev.off()
cat("\nMotifs are successfully plotted\n\n\n")

