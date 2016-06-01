#Want to plot the median copy number of each target for panel and wgs. 
#For wgs I will include some Mbp extra on each side to get rid of the noise. 
setwd("/home/Crisp/rebecka/")

#Find which bins belongs to which targets
targetBED <- read.table("/home/Crisp/rebecka/breastv7_QDNAseq_reb/BED_files/klevebring_clinseq_v3.targets.slopped.bed",sep="\t")
colnames(targetBED)<-c("chromosome","start","end","V4","V5","V6")
files_panel <- dir(pattern="panel.*qdnaseq.txt", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", recursive = TRUE)
files_wgs <- dir(pattern="wgs.*qdnaseq.txt", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", recursive = TRUE)
paneldata <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", files_panel[1], sep=""), header = TRUE, as.is=TRUE)
binsAtTarget <- paneldata[paneldata$use==TRUE, ] # these are the bins overlapping with targets
binsAtTarget$target <- NA

#targetBEDnoXYMT <- targetBED[-which(targetBED$chromosome=="MT"),]
#targetBEDnoXYMT <- targetBEDnoXYMT[-which(targetBEDnoXYMT$chromosome=="Y"),]
#targetBEDnoXYMT <- targetBEDnoXYMT[-which(targetBEDnoXYMT$chromosome=="X"),]

for (i in 1:length(binsAtTarget$start)) {
  for (j in 1:length(targetBED$start)) {
    if (!is.na(binsAtTarget$target[i])) {break}
    if (j==length(targetBED$start)) {print("didn't break")}
    if (binsAtTarget$chromosome[i]==targetBED$chromosome[j] && binsAtTarget$start[i]<=targetBED$start[j] && binsAtTarget$start[i+1]>targetBED$start[j]) {
      binsAtTarget$target[i] <- j 
    }
    if (binsAtTarget$chromosome[i]==targetBED$chromosome[j] && binsAtTarget$start[i]>targetBED$start[j] && binsAtTarget$start[i]<targetBED$end[j]) {
      binsAtTarget$target[i] <- j
    }
    }
    
  }
}

binsAtTarget_incomplete <- binsAtTarget
binsAtTarget_incomplete2 <- binsAtTarget

#Calculate the median copy number of each target


