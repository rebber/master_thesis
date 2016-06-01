setwd("/home/Crisp/rebecka/") #jag kan ej skriva i mappen /home/Crsip/rebecka/breastv7_QDNAseq_reb/ härifrån eftersom jag här är inne på vagrant användaren, får skriva i rebecka-mappen och sedan flytta till denna mapp från terminalen med min användare
biocLite('QDNAseq.hg19')
library(QDNAseq)
library(tools)
bins <- getBinAnnotations(binSize = 15)
bins
files_panel <- dir(pattern="T.*panel.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_panel <- grep(pattern="bai", files_panel, invert=TRUE, value=TRUE)
dir_panel <- dirname(files_panel)
files_wgs <- dir(pattern="T.*wgs.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_wgs <- grep(pattern="bai", files_wgs, invert=TRUE, value=TRUE)
dir_wgs <- dirname(files_wgs)

#QDNAseq for wgs
mybamfiles=paste("/home/Crisp/BREASTv7_bams_for_Rebecka/",files_wgs,sep="")
readCounts <- binReadCounts(bins,bamfiles=mybamfiles)
plot(readCounts,logTransform=FALSE,ylim=c(-50,200))
highlightFilters(readCounts,logTransform=FALSE,residual = TRUE, blacklist = TRUE)
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
pdf("isobarplot_wgs.pdf")
isobarPlot(readCountsFiltered)
dev.off()
readCountsFiltered <-estimateCorrection(readCountsFiltered)
pdf("noiseplot_wgs.pdf")
noisePlot(readCountsFiltered)
#dev.copy(pdf, "noiseplot_wgs.pdf")
dev.off()
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth,transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
plot(copyNumbersSegmented)
#tune segmentation?
copyNumbersCalled <- callBins(copyNumbersSegmented)
#filenames <- paste(basename(file_path_sans_ext(files_wgs)),".pdf",sep="")
pdf("lossgainplot_wgs.pdf")
plot(copyNumbersCalled)
dev.off()
cgh <- makeCgh(copyNumbersCalled)
cgh
#Här är QDNAseq klart

#Vill nu spara all data som txt för varje sample för sig
featureData <- copyNumbersCalled@featureData@data

for (i in 1:10) {
  sample <- featureData
  sample$readcount <- readCounts@assayData$counts[ ,i]
  sample$copynumber <- copyNumbersCalled@assayData$copynumber[ ,i]
  sample$segmented <- copyNumbersCalled@assayData$segmented[ ,i]
  sample$call <- copyNumbersCalled@assayData$calls[ ,i]
  sample$probdloss <- copyNumbersCalled@assayData$probdloss[ ,i]
  sample$probloss <- copyNumbersCalled@assayData$probloss[ ,i]
  sample$probnorm <- copyNumbersCalled@assayData$probnorm[ ,i]
  sample$probgain <- copyNumbersCalled@assayData$probgain[ ,i]
  sample$probamp <- copyNumbersCalled@assayData$probamp[ ,i]
  #head(sample)
  name <- paste(copyNumbersCalled@phenoData@data$name[i], "_reb.qdnaseq.txt",sep="")
  write.table(sample, name, col.names=TRUE, dec=".", quote=FALSE, sep="\t", row.names=FALSE)
}

#############################################################

#QDNAseq for panel
mybamfiles=paste("/home/Crisp/BREASTv7_bams_for_Rebecka/",files_panel,sep="")
readCounts <- binReadCounts(bins,bamfiles=mybamfiles)
plot(readCounts,logTransform=FALSE,ylim=c(-50,200))
highlightFilters(readCounts,logTransform=FALSE,residual = TRUE, blacklist = TRUE)
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE)
pdf("isobarplot_panel.pdf")
isobarPlot(readCountsFiltered)
dev.off()
readCountsFiltered <-estimateCorrection(readCountsFiltered)
pdf("noiseplot_panel.pdf")
noisePlot(readCountsFiltered)
#dev.copy(pdf, "noiseplot_panel.pdf")
dev.off()
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
plot(copyNumbersSmooth)
copyNumbersSegmented <- segmentBins(copyNumbersSmooth,transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
pdf("segmented_panel.pdf")
plot(copyNumbersSegmented)
dev.off()
#tune segmentation?
copyNumbersCalled <- callBins(copyNumbersSegmented)
#filenames <- paste(basename(file_path_sans_ext(files_panel)),".pdf",sep="")
pdf("lossgainplot_panel.pdf",fg="black")
plot(copyNumbersCalled)
dev.off()
cgh <- makeCgh(copyNumbersCalled)
cgh
#Här är QDNAseq klart

#Vill nu spara all data som txt för varje sample för sig
featureData <- copyNumbersCalled@featureData@data

for (i in 1:10) {
  sample <- featureData
  sample$readcount <- readCounts@assayData$counts[ ,i]
  sample$copynumber <- copyNumbersCalled@assayData$copynumber[ ,i]
  sample$segmented <- copyNumbersCalled@assayData$segmented[ ,i]
  sample$call <- copyNumbersCalled@assayData$calls[ ,i]
  sample$probdloss <- copyNumbersCalled@assayData$probdloss[ ,i]
  sample$probloss <- copyNumbersCalled@assayData$probloss[ ,i]
  sample$probnorm <- copyNumbersCalled@assayData$probnorm[ ,i]
  sample$probgain <- copyNumbersCalled@assayData$probgain[ ,i]
  sample$probamp <- copyNumbersCalled@assayData$probamp[ ,i]
  #head(sample)
  name <- paste(copyNumbersCalled@phenoData@data$name[i], "_reb.qdnaseq.txt",sep="")
  write.table(sample, name, col.names=TRUE, dec=".", quote=FALSE, sep="\t", row.names=FALSE)
}


# cgh1<-cgh[ ,1]
# exportBins(copyNumbersCalled, file="breastv7_QDNAseq_reb_selected.txt")
# head(copyNumbersCalled)
# exportBins(cgh, file="breastv7_QDNAseq_reb_selected2.txt",type=c("copynumber","segments","calls"))
# exportBins(cgh[ ,1], file="breastv7_QDNAseq_reb_selected_sample1.txt")
# dat <- read.table("breastv7_QDNAseq_reb_selected_sample1.txt", header = TRUE, as.is=TRUE)
# head(dat)



