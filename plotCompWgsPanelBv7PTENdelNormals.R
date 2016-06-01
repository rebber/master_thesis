#Plot comparison of the segments and bins for normals from CNVkit for panel data and for tumors from QDNAseq for wgs data, for Bv7 samples with PTEN homozygous deletion
#created by Rebecka Bergström on 2016-02-26
library(tools)
setwd("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl//")

files_cnr <- dir(pattern="cnr", "/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/batchResults/", recursive = TRUE)
files_cns <- dir(pattern="cns", "/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/batchResults/", recursive = TRUE)
files_wgs <- c("PTENdel1/PTENdel1T/wgs/PTENdel1N_wgs.qdnaseq.txt","PTENdel2/PTENdel2T/wgs/PTENdel2N_wgs.qdnaseq.txt","PTENdel3/PTENdel3T/wgs/PTENdel3N_wgs.qdnaseq.txt","PTENdel4/PTENdel4T/wgs/PTENdel4N_wgs.qdnaseq.txt")
#exchange PTENdel1 etc towards actual sample IDs

#Find the start and end position of PTEN on chr 10
genes_bed <- read.table("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/example_files.qdnaseq.genes.bed",as.is=TRUE)
PTEN_ind <- grep("PTEN",genes_bed$V4) #Finds both PTEN and PTENP1
PTEN_ind <- PTEN_ind[-which(PTEN_ind==grep("PTENP1",genes_bed$V4))] 
PTEN_bed <- genes_bed[PTEN_ind,] #The bed-data for PTEN
colnames(PTEN_bed) <- c("chromosome", "start", "end","gene","CN","strand")
PTEN_genestart <- PTEN_bed$start
PTEN_geneend <- PTEN_bed$end

#initiate lists for significance test results and bins used for significance testing
test1_panel_collection <- list()
test2_panel_collection <- list()
test1_wgs_collection <- list()
test2_wgs_collection <- list()
PTENbins_panel_collection <- list()
surrPTEN1_panel_collection <- list()
surrPTEN2_panel_collection <- list()
PTENbins_wgs_collection <- list()
surrPTEN1_wgs_collection <- list()
surrPTEN2_wgs_collection <- list()

#loop for significance testing and plotting of each sample
for (i in 1:length(files_cnr)) {
  dat_cnr <- read.table(paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/batchResults/", files_cnr[i], sep=""), header = TRUE, as.is=TRUE)
  chr10_cnr <- dat_cnr[which(dat_cnr$chromosome==10),] # choose the bins in chr 10
  dat_cns <- read.table(paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/batchResults/", files_cns[i], sep=""), header = TRUE, as.is=TRUE)
  chr10_cns <- dat_cns[which(dat_cns$chromosome==10),] # choose the segments in chr 10
  dat_wgs <- read.table(paste("/home/Crisp/clinseq/BREASTv7/",files_wgs[i],sep=""),header=TRUE,as.is=TRUE)
  chr10_wgs <- dat_wgs[which(dat_wgs$chromosome==10),]
  
  #collect data for PTEN and surrondings (5 Mb at each side, centered around 85 mb and 95 Mb, for use in significance test
  PTEN_startBin_panel <- which(chr10_cnr$start<=PTEN_genestart & chr10_cnr$end>=PTEN_genestart) #the index of the 1st bin overlapping with PTEN
  PTEN_endBin_panel <- which(chr10_cnr$start<=PTEN_geneend & chr10_cnr$end>=PTEN_geneend) #the index of the last bin overlapping with PTEN
  PTENbins_panel <- chr10_cnr[PTEN_startBin_panel:PTEN_endBin_panel,] #cnr-data for the bins overlapping with PTEN
  surrPTENstart1_panel <- which(chr10_cnr$start<=(82.5e6) & chr10_cnr$end>=(82.5e6)) #index of the 1st bins to use as control (aka suronding -> surr) 1 (left) in test
  surrPTENend1_panel <- which(chr10_cnr$start<=(87.5e6) & chr10_cnr$end>=(87.5e6)) #index of last bin to use as control 1 in test
  surrPTEN1_panel <- chr10_cnr[surrPTENstart1_panel:surrPTENend1_panel,] #cnr-data for left ctrl-region
  surrPTENstart2_panel <- which(chr10_cnr$start<=(92.5e6) & chr10_cnr$end>=(92.5e6)) #index of the 1st bins to use as control 2 (right) in test
  surrPTENend2_panel <- which(chr10_cnr$start<=(97.5e6) & chr10_cnr$end>=(97.5e6)) #index of last bin to use as control 2 in test
  surrPTEN2_panel <- chr10_cnr[surrPTENstart2_panel:surrPTENend2_panel,] #cnr-data for right ctrl-region  
  PTEN_startBin_wgs <- ceiling(PTEN_genestart/15000) #the index of the 1st bin overlapping with PTEN
  PTEN_endBin_wgs <- ceiling(PTEN_geneend/15000) #the index of the last bin overlapping with PTEN
  PTENbins_wgs <- chr10_wgs[PTEN_startBin_wgs:PTEN_endBin_wgs,] #wgs-data for the bins overlapping with PTEN
  surrPTEN1_wgs <- chr10_wgs[round(82.5e6/15e3):round(87.5e6/15e3),] #the left control region for wgs, 5 Mb centered around 85 Mb
  surrPTEN2_wgs <- chr10_wgs[round(92.5e6/15e3):round(97.5e6/15e3),] #the right control region for wgs, 5 Mb centered around 95 Mb
  
  #perform significance test
  test1_panel <- wilcox.test(PTENbins_panel$log2,surrPTEN1_panel$log2,alternative="l")
  test2_panel <- wilcox.test(PTENbins_panel$log2,surrPTEN2_panel$log2,alternative="l")
  test1_wgs <- wilcox.test(log2(PTENbins_wgs$copynumber),log2(surrPTEN1_wgs$copynumber),alternative="l")
  test2_wgs <- wilcox.test(log2(PTENbins_wgs$copynumber),log2(surrPTEN2_wgs$copynumber),alternative="l")
  
  #store the results and the bins used for testing
  test1_panel_collection[[i]] <- test1_panel
  test2_panel_collection[[i]] <- test2_panel
  test1_wgs_collection[[i]] <- test1_wgs
  test2_wgs_collection[[i]] <- test2_wgs
  PTENbins_panel_collection[[i]] <- PTENbins_panel
  surrPTEN1_panel_collection[[i]] <- surrPTEN1_panel
  surrPTEN2_panel_collection[[i]] <- surrPTEN2_panel
  PTENbins_wgs_collection[[i]] <- PTENbins_wgs
  surrPTEN1_wgs_collection[[i]] <- surrPTEN1_wgs
  surrPTEN2_wgs_collection[[i]] <- surrPTEN2_wgs
  
  
  #the indeces and coordinates for zoom in on PTEN 
  #for wgs
  startZoomIn <- floor(PTEN_genestart/15000)-50 #index of the first bin included in zoom in
  endZoomIn <- floor(PTEN_geneend/15000)+50 #index of last bin included in zoom in
  #for the panel
  zoomIn_cnr <- which(chr10_cnr$start>chr10_wgs$start[startZoomIn] & chr10_cnr$start<chr10_wgs$start[endZoomIn]) #which bins start in the zoom in window
  zoomIn_cns <- which(chr10_cns$start>chr10_wgs$start[startZoomIn] & chr10_cns$start<chr10_wgs$start[endZoomIn]) #which segments start in the zoom in window
  segmStartXcoordZoomIn <- c(chr10_wgs$start[startZoomIn],chr10_cns$start[zoomIn_cns]) #x-coordinates for segment starts in zoom in 
  segmEndXcoordZoomIn <- c(chr10_cns$end[zoomIn_cns-1],chr10_wgs$end[endZoomIn]) #x-coordinates for segment ends in zoom in
  if (length(zoomIn_cns)>0) {
    segmYcoordZoomIn <- chr10_cns$log2[c((zoomIn_cns[1]-1),zoomIn_cns)] #Y-coordinates for the segments in the zoom in
  } else {
    segmYcoordZoomIn <- chr10_cns$log2[which(chr10_cns$start<=chr10_wgs$start[startZoomIn] & chr10_cns$end>=chr10_wgs$start[startZoomIn])]
  }
  
  title_wgs <- paste("BREASTv7 log2 copy number chr 10 NORMALS ",basename(file_path_sans_ext(files_wgs[i])),sep="")
  title_wgs_zoomIn <- "zoom in on PTEN, wgs"
  title_panel <- paste("BREASTv7 log2 copy number chr 10 NORMALS ",basename(file_path_sans_ext(files_cnr[i])),".cnvkit",sep="")
  title_panel_zoomIn <- "zoom in on PTEN, panel"
  filename=paste("Bv7_PTENdel_compWgsPanel_Normals_",basename(file_path_sans_ext(file_path_sans_ext(gsub("_wgs","",files_wgs[i])))),".pdf",sep="")
  
  pdf(filename,width=20,height=9.8)
  layout(matrix(c(1,1,2,3,3,4), 2, 3, byrow = TRUE)) #plot for wgs and panel in the same graph
  
  #plot the wgs data from QDNAseq
  plot(chr10_wgs$start, log2(chr10_wgs$copynumber), xaxt="n",pch=16, col='#00000090', cex=0.3, ylim=c((-3.5),3.5), main=title_wgs,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN wgs")
  points(chr10_wgs$start, log2(chr10_wgs$segmented), cex=0.6, pch=16, col="#66CD00")
  abline(v=c(PTEN_genestart), h=0, lwd=0.5) #vertical line for start of PTEN and horizontal line for neutral CN
  abline(v=c(chr10_wgs$start[startZoomIn],chr10_wgs$start[endZoomIn]),lty=3)
  axis(1, at=seq(0,1.4e+08,20000000))
  text(0,3,paste("p-value PTEN vs left ctrl-reg (5 Mb around 85 Mb): ",signif(test1_wgs$p.value,digits=3),sep=""),pos=4)
  text(0,2.5,paste("p-value PTEN vs right ctrl-reg (5 Mb around 95 Mb): ",signif(test2_wgs$p.value,digits=3),sep=""),pos=4)
  #zoom in on PTEN
  plot(chr10_wgs$start[(startZoomIn):(endZoomIn)],log2(chr10_wgs$copynumber[(startZoomIn):(endZoomIn)]), xaxt="n",pch=16, col='#000000', cex=0.3, xlim=(c(chr10_wgs$start[(startZoomIn)]-1e4,chr10_wgs$start[(endZoomIn)]+1e4)), ylim=c((-3.5),3.5), main=title_wgs_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN wgs")
  segments(chr10_wgs$start[(startZoomIn):(endZoomIn)],log2(chr10_wgs$segmented[(startZoomIn):(endZoomIn)]),chr10_wgs$end[(startZoomIn):(endZoomIn)],log2(chr10_wgs$segmented[(startZoomIn):(endZoomIn)]), lwd=4, col="#66CD00")
  abline(v=c(PTEN_genestart,PTEN_geneend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=seq(8.9e7,9.05e7,0.05e7))
  
  #plot the panel data from CNVkit
  plot(chr10_cnr$start, chr10_cnr$log2, xaxt="n", pch=16, col='#00000099', cex=0.3, ylim=c((-3.5),3.5), main=title_panel,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN panel")
  segments(chr10_cns$start, chr10_cns$log2, chr10_cns$end, chr10_cns$log2, lwd=4, col="#66CD00")
  abline(v=PTEN_genestart, h=0, lwd=0.5) #vertical lines for the chromosomes and horizontal for neutral log2CN
  abline(v=c(chr10_wgs$start[startZoomIn],chr10_wgs$start[endZoomIn]),lty=3)
  axis(1, at=seq(0,1.4e+08,20000000))
  text(0,3,paste("p-value PTEN vs left ctrl-reg (5 Mb around 85 Mb): ",signif(test1_panel$p.value,digits=3),sep=""),pos=4)
  text(0,2.5,paste("p-value PTEN vs right ctrl-reg (5 Mb around 95 Mb): ",signif(test2_panel$p.value,digits=3),sep=""),pos=4)
  #zoom in on PTEN 
  plot(chr10_cnr$start[zoomIn_cnr], chr10_cnr$log2[zoomIn_cnr], xaxt="n", pch=16, col='#000000', cex=0.3, xlim=(c(chr10_wgs$start[(startZoomIn)]-1e4,chr10_wgs$start[(endZoomIn)]+1e4)), ylim=c((-3.5),3.5), main=title_panel_zoomIn,xlab="Position on Chr 10",cex.main=0.9, ylab="log2 CN panel")
  segments(segmStartXcoordZoomIn,segmYcoordZoomIn,segmEndXcoordZoomIn,segmYcoordZoomIn, lwd=4, col="#66CD00")
  abline(v=c(PTEN_genestart,PTEN_geneend), h=0, lwd=0.5) #vertical lines for start and end of PTEN and horizontal line for neutral CN
  axis(1,at=seq(8.9e7,9.05e7,0.05e7))
  
  dev.off()
}

#save the results from sign testing and the bins used for that 
saveRDS(test1_panel_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","test1_panelN.rds",sep=""))
saveRDS(test2_panel_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","test2_panelN.rds",sep=""))
saveRDS(PTENbins_panel_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","PTENbins_panelN.rds",sep=""))
saveRDS(surrPTEN1_panel_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","surrPTEN1_panelN.rds",sep=""))
saveRDS(surrPTEN2_panel_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","surrPTEN2_panelN.rds",sep=""))

saveRDS(test1_wgs_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","test1_wgsN",".rds",sep=""))
saveRDS(test2_wgs_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","test2_wgsN",".rds",sep=""))
saveRDS(PTENbins_wgs_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","PTENbins_wgsN",".rds",sep=""))
saveRDS(surrPTEN1_wgs_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","surrPTEN1_wgsN",".rds",sep=""))
saveRDS(surrPTEN2_wgs_collection,paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/","surrPTEN2_wgsN",".rds",sep=""))
