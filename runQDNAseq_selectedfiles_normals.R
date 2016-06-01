#QDNAseq for normals with panel


setwd("/home/Crisp/rebecka/") #jag kan ej skriva i mappen /home/Crsip/rebecka/breastv7_QDNAseq_reb/ härifrån eftersom jag här är inne på vagrant användaren, får skriva i rebecka-mappen och sedan flytta till denna mapp från terminalen med min användare
source("http://bioconductor.org/biocLite.R")
biocLite('QDNAseq.hg19')
biocLite('QDNAseq.hg19.1kbp.PE100')
library(QDNAseq)
library(tools)

# Get bins
#bins <- getBinAnnotations(binSize = 1)
#bedfiles <- c("clinseq_v3_targets_inverted_6col.bed","dukeExcludeRegions.bed","consensusBlacklist.bed")
#bedfiles <- c("dukeExcludeRegions.bed","consensusBlacklist.bed") #default
#bedfiles <- paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/BED_files/",bedfiles,sep="")
#bins@data$blacklist <- calculateBlacklist(bins@data, bedFiles=bedfiles) #använd för att läsa in blacklistning
#saveRDS(bins, "bins_1k_sr50_nontargbl.rds") #save the adjusted blacklist
#bins1ksr50ntbl<-readRDS("bins_1k_sr50_nontargbl.rds") # bins för panelen (om PE100 lyckas (nästa rad) ta dem istället)
#bins1kpe100 <- getBinAnnotations(binSize=1, type="PE100") #doesn't work, gets "Error in gzfile(file, "rb") : cannot open the connection" etc, verkar inte finnas nån fil för PE100
#bins1ksr50normbl <- getBinAnnotations(binSize = 1) # bins för wgs
#bins1ksr50_noblacklist <- bins1ksr50ntbl
#bins1ksr50_noblacklist@data$blacklist <- NULL #bins utan blacklist
#binstargetoverlap <- calculateBlacklist(bins1ksr50_noblacklist@data, bedFiles="/home/Crisp/rebecka/breastv7_QDNAseq_reb/BED_files/klevebring_clinseq_v3.targets.slopped.bed") #beräkna varje bins överlapp med targetregioner
#saveRDS(binstargetoverlap, "1kbpbinstargetoverlap.rds") #flyttad till breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/
binstargetoverlap <- readRDS("breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/1kbpbinstargetoverlap.rds")
binsAtTarget <- which(binstargetoverlap>0) #indeces for which bins use=TRUE should be set (we want to use all bins overlapping with target)
rm(binstargetoverlap)


#Get bam files
files_panel <- dir(pattern="N.*panel.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_panel <- grep(pattern="bai", files_panel, invert=TRUE, value=TRUE)
files_wgs <- dir(pattern="N.*wgs.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_wgs <- grep(pattern="bai", files_wgs, invert=TRUE, value=TRUE)




#############################################################

#QDNAseq for panel
#bins1ksr50ntbl<-readRDS("breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/bins_1k_sr50_nontargbl.rds") # bins för panelen (om PE100 lyckas ta dem istället)
mybamfiles=paste("/home/Crisp/BREASTv7_bams_for_Rebecka/",files_panel,sep="")
readCounts <- binReadCounts(bins1ksr50ntbl,bamfiles=mybamfiles)
saveRDS(readCounts,"readCounts_panel_normal.rds") 
rm(bins1ksr50ntbl)
readCounts <- readRDS("breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/readCounts_panel_normal.rds")
##################################################
#Current filtration: Use all bins overlapping w/ target, no others
#################################################
readCountsFiltered <- applyFilters(readCounts, residual = TRUE, blacklist = TRUE) #important since it removes X and Y chromosomes, the blacklisting is not important since the use is changed manually in next steps
#Enligt Markus borde inte de default svartlistade regionerna överlappa med våra targets
readCountsFiltered@featureData@data$use <- FALSE
readCountsFiltered@featureData@data$use[binsAtTarget] <- TRUE
readCountsFiltered <-estimateCorrection(readCountsFiltered)
pdf("noiseplot_panel_normal.pdf")
noisePlot(readCountsFiltered)
dev.off()
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
#rm(readCountsFiltered,copyNumbers,copyNumbersNormalized)
#copyNumbersSegmented <- segmentBins(copyNumbersSmooth,transformFun="sqrt")
#copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
#tune segmentation?
#copyNumbersCalled <- callBins(copyNumbersSegmented)
#cgh <- makeCgh(copyNumbersCalled)
#cgh
#Här är QDNAseq klart

#Vill nu spara all data som txt för varje sample för sig
featureData <- copyNumbersSmooth@featureData@data

for (i in 1:10) {
  sample <- featureData
  sample$readcount <- readCounts@assayData$counts[ ,i]
  sample$copynumber <- copyNumbersSmooth@assayData$copynumber[ ,i]
  #   sample$segmented <- copyNumbersCalled@assayData$segmented[ ,i]
  #   sample$call <- copyNumbersCalled@assayData$calls[ ,i]
  #  sample$probdloss <- copyNumbersCalled@assayData$probdloss[ ,i]
  #  sample$probloss <- copyNumbersCalled@assayData$probloss[ ,i]
  #  sample$probnorm <- copyNumbersCalled@assayData$probnorm[ ,i]
  #  sample$probgain <- copyNumbersCalled@assayData$probgain[ ,i]
  #  sample$probamp <- copyNumbersCalled@assayData$probamp[ ,i]
  #head(sample)
  name <- paste(copyNumbersSmooth@phenoData@data$name[i], "_rebadj.qdnaseq.txt",sep="")
  write.table(sample, name, col.names=TRUE, dec=".", quote=FALSE, sep="\t", row.names=FALSE)
}

#remove the QDNAseq objects for panel
"rm(readCounts,copyNumbersCalled,copyNumbersSegmented,copyNumbersSmooth,featureData,sample)
