#Match log2 copy number ratios for genes from wgs with QDNAseq and panel with CNVkit. Plot scatter plots. 
#Here for MSI-samples
library(MASS)
library(tools)
setwd("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/Rplots/scatter_panelVSwgs")

files_wgs <- dir(pattern="genes.bed", "/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/", recursive = TRUE)
files_panel <- dir(pattern="genes.bed","/home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/", recursive=TRUE)

counter_list <- c()
log2CN_panel_all <- c()
log2CN_wgs_all <- c()
log2CN_panel_remextr_all <- c()
log2CN_wgs_remextr_all <- c()
fits<-rep(NA,length(files_wgs))
fits_remextr<-rep(NA,length(files_wgs))
for (i in 1:length(files_wgs)) {
  data_wgs <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/QDNAseqWGS/segments2genes/",files_wgs[i],sep=""),as.is=TRUE)
  colnames(data_wgs) <- c("chromosome", "start", "end", "gene", "CN", "strand")
  data_wgs$log2 <- log2(data_wgs$CN)
  #data_wgs <- data_wgs[-which(is.na(data_wgs$log2)),]
  data_panel <- read.table(paste("/home/rebecka.bergstrom/cnvkit/MSI-PILOT/segments2genes/",files_panel[i],sep=""), as.is=TRUE)
  colnames(data_panel) <- c("chromosome", "start", "end", "gene", "log2", "strand")
  
  #remove extreme values (log2CN<=-3.5 and log2CN>=3.5)
  extrrem <- which(data_panel$log2<=(-3.5))
  extrrem <- c(extrrem, which(data_wgs$log2<=(-3.5)))
  extrrem <- c(extrrem, which(data_wgs$log2>=(3.5)))
  extrrem <- c(extrrem, which(data_panel$log2>=(3.5)))
  if (length(extrrem>0)) {
    log2CN_panel_remextr <- data_panel$log2[-extrrem]
    log2CN_wgs_remextr <- data_wgs$log2[-extrrem]
  } else {
    log2CN_panel_remextr <- data_panel$log2
    log2CN_wgs_remextr <- data_wgs$log2
  }
  
  #fit a line both with and without extreme values removed
  fit <- rlm(data_panel$log2 ~ data_wgs$log2)
  fit_remextr <- rlm(log2CN_panel_remextr ~ log2CN_wgs_remextr)
  
  #plot
  title1 <- paste("log2CN panel vs log2CN wgs per gene ",basename(file_path_sans_ext(files_panel[i])),sep="")
  title2 <- paste("log2CN panel vs log2CN wgs per gene without extreme points ",basename(file_path_sans_ext(files_panel[i])),sep="")
  filename <- paste("scatter_", basename(file_path_sans_ext(files_panel[i])), ".pdf",sep="")
  
  pdf(filename,width=10,height=10)
  plot(data_wgs$log2,data_panel$log2,pch=16, col='#00000050',cex=0.5, main=title1, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
  abline(coef(fit))
  abline(c(fit$coefficients[1],1), col=4)
  legend("topleft",legend=c("data", paste("fitted line k=",round(fit$coefficients[2],digits=2)), "k=1"),col=c('#00000050',153,4),pch=c(16,NA,NA),lwd=c(NA,1,1))
  plot(log2CN_wgs_remextr,log2CN_panel_remextr, pch=16, col='#00000050',cex=0.5, main=title2, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
  abline(coef(fit_remextr))
  abline(c(fit_remextr$coefficients[1],1), col=4)
  legend("topleft",legend=c("data", paste("fitted line k=",round(fit_remextr$coefficients[2],digits=2)), "k=1"),col=c('#00000050',153,4),pch=c(16,NA,NA),lwd=c(NA,1,1))
  dev.off()
  
  #How many genes matched (i.e. had a log2CN value both for panel and wgs)
  counter <- sum(!is.na(data_wgs$log2) & !is.na(data_panel$log2))
  print(paste("Matched",counter,"genes"))
  
  #save entries from each sample to analyze outside for-loop
  counter_list[[i]] <- counter
  fits[[i]] <- fit$coefficients[2]
  fits_remextr[[i]] <- fit_remextr$coefficients[2]
  log2CN_panel_all <- c(log2CN_panel_all, data_panel$log2)
  log2CN_wgs_all <- c(log2CN_wgs_all, data_wgs$log2)
  log2CN_panel_remextr_all <- c(log2CN_panel_remextr_all, log2CN_panel_remextr)
  log2CN_wgs_remextr_all <- c(log2CN_wgs_remextr_all, log2CN_wgs_remextr)
  
}


#find the fit and plot the concatenation from all samples
fit_all <- rlm(log2CN_panel_all ~ log2CN_wgs_all)
fit_remextr_all <- rlm(log2CN_panel_remextr_all ~ log2CN_wgs_remextr_all)
title1 <- "log2CN panel vs log2CN wgs per gene, all MSI samples together"
title2 <- "log2CN panel vs log2CN wgs per gene without extreme points, all MSI samples together"
filename <- "scatter_allMSIsamples.genes.pdf"

pdf(filename,width=10,height=10)
plot(log2CN_wgs_all,log2CN_panel_all,pch=16, col='#00000050',cex=0.5, main=title1, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
abline(coef(fit_all))
abline(c(fit_all$coefficients[1],1), col=4)
legend("topleft",legend=c("data", paste("fitted line k=",round(fit_all$coefficients[2],digits=2)), "k=1"),col=c('#00000050',153,4),pch=c(16,NA,NA),lwd=c(NA,1,1))
plot(log2CN_wgs_remextr_all,log2CN_panel_remextr_all, pch=16, col='#00000050',cex=0.5, main=title2, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
abline(coef(fit_remextr_all))
abline(c(fit_remextr_all$coefficients[1],1), col=4)
legend("topleft",legend=c("data", paste("fitted line k=",round(fit_remextr_all$coefficients[2],digits=2)), "k=1"),col=c('#00000050',153,4),pch=c(16,NA,NA),lwd=c(NA,1,1))
dev.off()



