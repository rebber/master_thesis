#Want to plot the median copy number of each target for panel and wgs. 
#For wgs I will include some Mbp extra on each side to get rid of the noise. 
setwd("/home/Crisp/rebecka/")
#install.packages("MASS")
library(MASS)
library(tools)

#Find which bins belongs to which targets
targetBED <- read.table("/home/Crisp/rebecka/breastv7_QDNAseq_reb/BED_files/klevebring_clinseq_v3.targets.slopped.bed",sep="\t")
colnames(targetBED)<-c("chromosome","start","end","V4","V5","V6")
files_panel_T <- dir(pattern="T.*panel.*qdnaseq.txt", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", recursive = TRUE)
files_wgs_T <- dir(pattern="T.*wgs.*qdnaseq.txt", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", recursive = TRUE)
paneldata <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", files_panel_T[1], sep=""), header = TRUE, as.is=TRUE)
binsAtTarget <- which(paneldata$use==TRUE) # these are the bin indeces overlapping with targets
fits<-rep(NA,10)

#Calculate median CN for each bin overlapping with targets and surronding bins for wgs
for (i in 1:length(files_wgs_T)) {
  dat_wgs <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", files_wgs_T[i], sep=""), header = TRUE, as.is=TRUE)
  datOverlap_wgs <- dat_wgs[binsAtTarget,]
  for (j in 1:length(binsAtTarget)) {
    datOverlap_wgs$medianCN[j] <- median(dat_wgs$copynumber[(binsAtTarget[j]-100):(binsAtTarget[j]+100)],na.rm=TRUE)
  }
  dat_panel <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_adj_alltargbins/", files_panel_T[i], sep=""), header = TRUE, as.is=TRUE)
  datOverlap_panel <- dat_panel[binsAtTarget,]
  log2CN_wgs <- log2(datOverlap_wgs$medianCN)
  NArem <- which(is.na(log2CN_wgs))
  log2CN_wgs <- log2CN_wgs[-NArem]
  log2CN_panel <- log2(datOverlap_panel$copynumber)
  log2CN_panel <- log2CN_panel[-NArem]
  lowrem <- which(log2CN_panel<=(-5))
  log2CN_panel_remlow <- log2CN_panel[-lowrem]
  log2CN_wgs_remlow <- log2CN_wgs[-lowrem]
  fit <- rlm(log2CN_panel_remlow ~ log2CN_wgs_remlow)
  title1 <- paste("log2CN panel vs median log2CN wgs ",basename(file_path_sans_ext(files_panel_T[i])),sep="")
  title2 <- paste("log2CN panel vs median log2CN wgs w/o low coverage ",basename(file_path_sans_ext(files_panel_T[i])),sep="")
  filename <- paste("scatter_", basename(file_path_sans_ext(files_panel_T[i])), ".pdf",sep="")
  pdf(filename,width=10,height=10)
  plot(log2CN_wgs, log2CN_panel, pch=16, col='#00000050',cex=0.5, main=title1, cex.main=0.8, xlim=c(-7,7), ylim=c(-7,7))
  plot(log2CN_wgs_remlow,log2CN_panel_remlow, pch=16, col='#00000050',cex=0.5, main=title2, cex.main=0.8,xlim=c(-5,6), ylim=c(-5,6))
  abline(coef(fit))
  abline(c(fit$coefficients[1],1), col=4)
  dev.off()
  
  fits[[i]] <- fit$coefficients[2] 
  
  #cor(log2CN_wgs, log2CN_panel)
  #fit <- lm(log2CN_panel ~ log2CN_wgs) 
  #fit <- rlm(log2CN_panel_remlow ~ log2CN_wgs_remlow)
  #fit2 <- line(log2CN_panel_remlow,log2CN_wgs_remlow)
  
}