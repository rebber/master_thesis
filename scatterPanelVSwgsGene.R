#Match log2 copy number ratios for genes from wgs with QDNAseq and panel with CNVkit. Plot scatter plots. 
library(MASS)
library(tools)
setwd("/home/rebecka.bergstrom/cnvkit/1stRunBreast/Rplots/scatter_panelVSwgs/")

files_wgs <- dir(pattern="genes.bed", "/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_nonadj/segments2genes", recursive = TRUE)
files_panel <- dir(patter="cns","/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun", recursive=TRUE)

counter_list <- list()
fits<-rep(NA,length(files_wgs))
fits_remlow<-rep(NA,length(files_wgs))
for (i in 1:length(files_wgs)) {
  data_wgs <- read.table(paste("/home/Crisp/rebecka/breastv7_QDNAseq_reb/QDNAseq_nonadj/segments2genes/",files_wgs[i],sep=""),as.is=TRUE)
  colnames(data_wgs) <- c("chromosome", "start", "end", "gene", "CN", "strand")
  data_wgs$log2 <- log2(data_wgs$CN)
  data_wgs <- data_wgs[-which(is.na(data_wgs$log2)),]
  data_panel <- read.table(paste("/home/rebecka.bergstrom/cnvkit/1stRunBreast/batchRun/",files_panel[i],sep=""), as.is=TRUE, header=TRUE)
  
  combined <- data_wgs
  combined$log2panel <- rep(NA,times=length(combined$log2))
  colnames(combined) <- c("chromosome", "start", "end", "gene", "CNpanel", "strand","log2wgs","log2panel")
  
  counter <- rep(NA,times=length(combined$log2))
  
  for (j in 1:length(data_wgs$gene)) {
    matching <- grep(data_wgs$gene[j],data_panel$gene)
    counter[j] <- length(matching)
    if (length(matching)>0) {
      log2panel <- data_panel$log2[matching]
      combined$log2panel[j] <- median(log2panel)
    }
  }
  print(paste("Matched",sum(counter>0),"genes"))
  counter_list[[i]] <- counter
  
  lowrem <- which(combined$log2panel<=(-3.5))
  log2CN_panel_remlow <- combined$log2panel[-lowrem]
  log2CN_wgs_remlow <- combined$log2wgs[-lowrem]
  fit <- rlm(combined$log2panel ~ combined$log2wgs)
  fit_remlow <- rlm(log2CN_panel_remlow ~ log2CN_wgs_remlow)
  fits[[i]] <- fit$coefficients[2]
  fits_remlow[[i]] <- fit_remlow$coefficients[2]
  title1 <- paste("log2CN panel vs log2CN wgs per gene ",basename(file_path_sans_ext(files_panel[i])),sep="")
  title2 <- paste("log2CN panel vs log2CN wgs per gene without low coverage points ",basename(file_path_sans_ext(files_panel[i])),sep="")
  filename <- paste("scatter_", basename(file_path_sans_ext(files_panel[i])), ".pdf",sep="")
  
  pdf(filename,width=10,height=10)
  plot(combined$log2wgs,combined$log2panel,pch=16, col='#00000050',cex=0.5, main=title1, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
  plot(log2CN_wgs_remlow,log2CN_panel_remlow, pch=16, col='#00000050',cex=0.5, main=title2, cex.main=0.8, xlab="log2 CN wgs", ylab="log2 CN panel")
  abline(coef(fit_remlow))
  abline(c(fit_remlow$coefficients[1],1), col=4)
  legend("topleft",legend=c("data", paste("fitted line k=",round(fit_remlow$coefficients[2],digits=2)), "k=1"),col=c('#00000050',153,4),pch=c(16,NA,NA),lwd=c(NA,1,1))
  dev.off()
  
  
}


counter_mat <- do.call(cbind,counter_list)
