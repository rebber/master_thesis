#make boxplots of log2CN for the bins in PTEN and the surronding control regions, for Bv7 samples with homozygous PTEN deletion
#created by Rebecka Bergström on 2016-02-26

library(tools)
setwd("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/boxplots/wilctest_5x2Mbctrl/")
samples <- c("PTENdel1","PTENdel2","PTENdel3","PTENdel4") #exchange for actual IDs of PTEN del samples
tumPerc_steps <- seq(100,0,-10)


#read the files containing the desired bins and their CN for pure tumors and normals, wgs- and panel-data
test1_panel100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test1_panelT.rds")
test2_panel100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test2_panelT.rds")
PTENbins_panel100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/PTENbins_panelT.rds")
surrPTEN1_panel100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN1_panelT.rds")
surrPTEN2_panel100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN2_panelT.rds")
test1_wgs100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test1_wgs.rds")
test2_wgs100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test2_wgs.rds")
PTENbins_wgs100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/PTENbins_wgs.rds")
surrPTEN1_wgs100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN1_wgs.rds")
surrPTEN2_wgs100 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/full_tumors/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN2_wgs.rds")
test1_panel0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test1_panelN.rds")
test2_panel0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test2_panelN.rds")
PTENbins_panel0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/PTENbins_panelN.rds")
surrPTEN1_panel0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN1_panelN.rds")
surrPTEN2_panel0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN2_panelN.rds")
test1_wgs0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test1_wgsN.rds")
test2_wgs0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/test2_wgsN.rds")
PTENbins_wgs0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/PTENbins_wgsN.rds")
surrPTEN1_wgs0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN1_wgsN.rds")
surrPTEN2_wgs0 <- readRDS("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/normals/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/surrPTEN2_wgsN.rds")

#read the files containing the desired bins and their CN for diluted panel data
for (i in 2:(length(tumPerc_steps)-1)) {
  tumPerc <- tumPerc_steps[i]
  path_rds <- paste("/home/rebecka.bergstrom/cnvkit/BREASTv7_PTEN_del/homozDelPTEN/panelV1_CORRECT/dilutions/tumorPerc_",tumPerc, "/Rplots/segmentPlotsComp/wilcoxontest/5x2Mbctrl/", sep="") #where the rds-files are to be found 
  assign(paste("test1_panel",tumPerc,sep=""),readRDS(paste(path_rds,"test1_panel",tumPerc,".rds",sep="")))
  assign(paste("test2_panel",tumPerc,sep=""),readRDS(paste(path_rds,"test2_panel",tumPerc,".rds",sep="")))
  assign(paste("PTENbins_panel",tumPerc,sep=""),readRDS(paste(path_rds,"PTENbins_panel",tumPerc,".rds",sep="")))
  assign(paste("surrPTEN1_panel",tumPerc,sep=""),readRDS(paste(path_rds,"surrPTEN1_panel",tumPerc,".rds",sep="")))
  assign(paste("surrPTEN2_panel",tumPerc,sep=""),readRDS(paste(path_rds,"surrPTEN2_panel",tumPerc,".rds",sep="")))
}

pvalues1_panel_all <- matrix(nrow=length(samples),ncol=length(tumPerc_steps))
pvalues2_panel_all <- matrix(nrow=length(samples),ncol=length(tumPerc_steps))
pvalues1_wgs_all <- matrix(nrow=length(samples),ncol=2)
pvalues2_wgs_all <- matrix(nrow=length(samples),ncol=2)

#create the boxplots sample wise
for (i in 1:length(samples)) {
  log2CNs_wgs <- list()
  log2CNs_panel <- list()
  pvalues1_wgs <- c()
  pvalues2_wgs <- c()
  pvalues1_panel <- c()
  pvalues2_panel <- c()
  #the bins and p-values for wgs (tumor & normal)
  log2CNs_wgs[[1]] <- log2(surrPTEN1_wgs100[[i]]$copynumber)
  log2CNs_wgs[[2]] <- log2(PTENbins_wgs100[[i]]$copynumber)
  log2CNs_wgs[[3]] <- log2(surrPTEN2_wgs100[[i]]$copynumber)
  log2CNs_wgs[[4]] <- log2(surrPTEN1_wgs0[[i]]$copynumber)
  log2CNs_wgs[[5]] <- log2(PTENbins_wgs0[[i]]$copynumber)
  log2CNs_wgs[[6]] <- log2(surrPTEN2_wgs0[[i]]$copynumber)
  pvalues1_wgs[1] <- test1_wgs100[[i]]$p.value
  pvalues2_wgs[1] <- test2_wgs100[[i]]$p.value
  pvalues1_wgs[2] <- test1_wgs0[[i]]$p.value
  pvalues2_wgs[2] <- test2_wgs0[[i]]$p.value
  #the bins and p-values for panel (from 100 % tumor via dilutions to 100 % normals)
  for (j in 1:length(tumPerc_steps)) {
    log2CNs_panel[[j*3-2]] <- eval(parse(text=paste("surrPTEN1_panel",tumPerc_steps[j],"[[i]]$log2",sep="")))
    log2CNs_panel[[j*3-1]] <- eval(parse(text=paste("PTENbins_panel",tumPerc_steps[j],"[[i]]$log2",sep="")))
    log2CNs_panel[[j*3]] <- eval(parse(text=paste("surrPTEN2_panel",tumPerc_steps[j],"[[i]]$log2",sep="")))
    pvalues1_panel[j] <- eval(parse(text=paste("test1_panel",tumPerc_steps[j],"[[i]]$p.value",sep="")))
    pvalues2_panel[j] <- eval(parse(text=paste("test2_panel",tumPerc_steps[j],"[[i]]$p.value",sep="")))
  }
  
  #store p values
  pvalues1_panel_all[i,] <- pvalues1_panel
  pvalues2_panel_all[i,] <- pvalues2_panel
  pvalues1_wgs_all[i,] <- pvalues1_wgs
  pvalues2_wgs_all[i,] <- pvalues2_wgs
  
  #plot
  title1 <-"panel data"
  title2 <- "WGS data"
  filename <- paste(samples[i],"_boxplot_dil_wilc5x2Mb.jpeg",sep="")
  #set at=xpos (see below) in boxplot() if boxew should be grouped together
  xpos <- c(1,1.8,2.6,4,4.8,5.6,7,7.8,8.6,10,10.8,11.6,13,13.8,14.6,16,16.8,17.6,19,19.8,20.6,22,22.8,23.6,25,25.8,26.6,28,28.8,29.6,31,31.8,32.6)+0.2 
  xpos_wgs <- c(1.2,2,2.8,4.2,5,5.8)
#  pdf(filename,width=20,height=9.8)
  jpeg(filename,width=9,height=6,quality=90,units="in",res=600)
  layout(matrix(c(1,1,1,2), 1, 4, byrow = TRUE))
  par(mar=c(4.7,4.3,3.8,1.7))
  boxplot(log2CNs_panel,col=c("dodgerblue4","gold","darkgreen"),ylim=c(-1.8,1.8),xaxt="n",ylab="log2 CN ratio",xlab="Part of original tumor content (%)",main=title1,at=xpos,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  abline(h=0,lty=2)
  abline(v=seq(0.5,(length(log2CNs_panel)+0.5),3),lty=3)
  legend("topright",c("left ctrl reg","PTEN","right ctrl reg"),fill=c("dodgerblue4","gold","darkgreen"))
  axis(1, at=seq(2,(length(log2CNs_panel)-1),3), labels=as.character(tumPerc_steps),cex.axis=1.5)
  text(0,1.2,"p-values PTEN vs left ctrl-reg:",pos=4,cex=1)
  text(seq(2,(length(log2CNs_panel)-1),3),rep(1.1,length(pvalues1_panel)/3),labels=signif(pvalues1_panel,digits=3),cex=1)
  text(0,0.9,"p-values PTEN vs right ctrl-reg:",pos=4,cex=1)
  text(seq(2,(length(log2CNs_panel)-1),3),rep(0.8,length(pvalues2_panel)/3),labels=signif(pvalues2_panel,digits=3),cex=1)
  boxplot(log2CNs_wgs,col=c("dodgerblue4","gold","darkgreen"),ylim=c(-1.8,1.8),xaxt="n",ylab="log2 CN ratio",xlab="Sample type",main=title2,at=xpos_wgs,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  abline(h=0,lty=2)
  abline(v=3.5, lty=3)
  legend("topright",c("left ctrl reg","PTEN","right ctrl reg"),fill=c("dodgerblue4","gold","darkgreen"))
  axis(1, at=c(2,5), labels=c("tumor","normal"),cex.axis=1.5)
#  text(0.5,1.7,"p-values PTEN vs left ctrl-reg:",pos=4,cex=1)
  text(c(2,5),c(1.1,1.1),labels=signif(pvalues1_wgs,3),cex=1)
#  text(0.5,1.2,"p-values PTEN vs right ctrl-reg:",pos=4,cex=1)
  text(c(2,5),c(0.8,0.8),labels=signif(pvalues2_wgs,3),cex=1)
  dev.off()

  #plot without wgs, also don't print p-values
  filename <- paste(samples[i],"_boxplot_dil_wilc5x2Mb_noPval_noWGS.jpeg",sep="")
  jpeg(filename,width=9,height=6,quality=90,units="in",res=600)
  layout(matrix(1,1,1, byrow = TRUE))
  par(mar=c(4.7,4.3,3.8,1.7))
  boxplot(log2CNs_panel,col=c("dodgerblue4","gold","darkgreen"),ylim=c(-1.8,1.8),xaxt="n",ylab="log2 CN ratio",xlab="Part of original tumor content (%)",main=title1,at=xpos,cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  abline(h=0,lty=2)
  abline(v=seq(0.5,(length(log2CNs_panel)+0.5),3),lty=3)
  legend("topright",c("left ctrl reg","PTEN","right ctrl reg"),fill=c("dodgerblue4","gold","darkgreen"))
  axis(1, at=seq(2,(length(log2CNs_panel)-1),3), labels=as.character(tumPerc_steps),cex.axis=1.5)
  #text(0,1.2,"p-values PTEN vs left ctrl-reg:",pos=4,cex=1)
  #text(seq(2,(length(log2CNs_panel)-1),3),rep(1.1,length(pvalues1_panel)/3),labels=signif(pvalues1_panel,digits=3),cex=1)
  #text(0,0.9,"p-values PTEN vs right ctrl-reg:",pos=4,cex=1)
  #text(seq(2,(length(log2CNs_panel)-1),3),rep(0.8,length(pvalues2_panel)/3),labels=signif(pvalues2_panel,digits=3),cex=1)
  dev.off()

}


#####################

#plot all p-values
filename_all_p_val <- "all_pVal.jpeg"
jpeg(filename_all_p_val,width=9,height=6,quality=90,units="in",res=600)
layout(matrix(1:2,1,2))
plot(tumPerc_steps,pvalues1_panel_all[1,],log="y",ylim=c(1e-10,1),type="n",xlab="Part of original tumor content (%)",ylab="p-value PTEN vs control regions", main="panel data",xlim=c(100,0))
points(tumPerc_steps,pvalues1_panel_all[1,],type="o",pch=8,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[1,],type="o",pch=8,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[2,],type="o",pch=15,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[2,],type="o",pch=15,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[3,],type="o",pch=16,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[3,],type="o",pch=16,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[4,],type="o",pch=17,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[4,],type="o",pch=17,lty=1,col="darkgreen")
abline(h=c(1e-3),lty=2)
legend_text <- c("sample:", 1:4, "ctrl reg:", "left","right", "cut-off")
legend("topleft",legend_text,pch=c(NA,8,15:17,rep(NA,4)),col=c(NA,rep(1,4),NA,"dodgerblue4","darkgreen",1),lty=c(rep(NA,6),1,1,2),lwd=c(rep(NA,6),2,2,1),cex=0.8)
plot(c(100,0),pvalues1_wgs_all[1,],log="y",ylim=c(1e-10,1),type="n",xlab="Part of original tumor content (%)",ylab="p-value PTEN vs control regions",main="WGS data",xlim=c(100,0))
points(c(100,0),pvalues1_wgs_all[1,],type="p",pch=8,col="dodgerblue4")
points(c(100,0),pvalues2_wgs_all[1,],type="p",pch=8,col="darkgreen")
points(c(100,0),pvalues1_wgs_all[2,],type="p",pch=15,col="dodgerblue4")
points(c(100,0),pvalues2_wgs_all[2,],type="p",pch=15,col="darkgreen")
points(c(100,0),pvalues1_wgs_all[3,],type="p",pch=16,col="dodgerblue4")
points(c(100,0),pvalues2_wgs_all[3,],type="p",pch=16,col="darkgreen")
points(c(100,0),pvalues1_wgs_all[4,],type="p",pch=17,col="dodgerblue4")
points(c(100,0),pvalues2_wgs_all[4,],type="p",pch=17,col="darkgreen")
#abline(h=c(1e-4,1e-3,1e-2),lty=2)
dev.off()

#no wgs
filename_all_p_val_noWGS <- "all_pVal_noWGS.jpeg"
jpeg(filename_all_p_val_noWGS,width=6,height=6,quality=90,units="in",res=600)
layout(matrix(1,1,1))
plot(tumPerc_steps,pvalues1_panel_all[1,],log="y",ylim=c(1e-10,1),type="n",xlab="Part of original tumor content (%)",ylab="p-value PTEN vs control regions",xlim=c(100,0))
points(tumPerc_steps,pvalues1_panel_all[1,],type="o",pch=8,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[1,],type="o",pch=8,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[2,],type="o",pch=15,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[2,],type="o",pch=15,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[3,],type="o",pch=16,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[3,],type="o",pch=16,lty=1,col="darkgreen")
points(tumPerc_steps,pvalues1_panel_all[4,],type="o",pch=17,lty=1,col="dodgerblue4")
points(tumPerc_steps,pvalues2_panel_all[4,],type="o",pch=17,lty=1,col="darkgreen")
abline(h=c(1e-3),lty=2)
legend_text <- c("sample:", 1:4, "ctrl reg:", "left","right", "cut-off")
legend("topleft",legend_text,pch=c(NA,8,15:17,rep(NA,4)),col=c(NA,rep(1,4),NA,"dodgerblue4","darkgreen",1),lty=c(rep(NA,6),1,1,2),lwd=c(rep(NA,6),2,2,1),cex=0.8)
dev.off()

pvalues1_panel_all[,11]/pvalues1_panel_all[,1]
pvalues1_wgs_all[,2]/pvalues1_wgs_all[,1]

pvalues2_panel_all[,11]/pvalues2_panel_all[,1]
pvalues2_wgs_all[,2]/pvalues2_wgs_all[,1]
