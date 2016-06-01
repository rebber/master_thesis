targets <- read.table("/home/Crisp/rebecka/breastv7_QDNAseq_reb/BED_files/klevebring_clinseq_v3.targets.slopped.bed")
targets <- targets[-which(targets[,1]=="MT"),]
targets <- targets[-which(targets[,1]=="Y"),]
targets <- targets[-which(targets[,1]=="X"),]
targetlengths <- targets[ ,3]-targets[ ,2]
barplot(sort(targetlengths))
boxplot(targetlengths)
quantile(targetlengths,probs=seq(0,1,0.1))
sum(targetlengths>1000)
hist(targets[ ,1])
head(targets[ ,1])
sum(targets[ ,1]=="MT")
smalltargets <- which(targetlengths<60) 
verysmalltargets <- which(targetlengths==40)
vsm_bedform <- targets[verysmalltargets,]
sm_bedform <- targets[smalltargets,]
write.table(vsm_bedform,file="targetlength40.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(sm_bedform,file="targetlength60.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
