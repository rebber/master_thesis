library(QDNAseq)
bins <- getBinAnnotations(binSize = 15)
bins
files_panel <- dir(pattern="T.*panel.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_panel <- grep(pattern="bai", files_panel, invert=TRUE, value=TRUE)
dir_panel <- dirname(files_panel)
files_wgs <- dir(pattern="T.*wgs.*.bam", "/home/Crisp/BREASTv7_bams_for_Rebecka", recursive = TRUE)
files_wgs <- grep(pattern="bai", files_wgs, invert=TRUE, value=TRUE)
dir_wgs <- dirname(files_wgs)

bamfile=paste("/home/Crisp/BREASTv7_bams_for_Rebecka/",files_wgs[1],sep="")
readCounts <- binReadCounts(bins,bamfiles=bamfile)


