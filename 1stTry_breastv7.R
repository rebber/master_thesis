# ---
#   title: "breastv7-wgs"
# author: "dakl"
# date: "26 January 2016"
# output: html_document
# ---
  
  # ClinSeq Breast
  
#   ```{r load-libs, echo=FALSE, warning=FALSE, results='hide', message=FALSE, cache=TRUE}
# ```
library(QDNAseq)
library(gplots)

files <- dir(pattern="qdnaseq.txt", "/home/Crisp/clinseq/BREASTv7/", recursive = TRUE)

#detta väljer ut fil #1, loopa och plotta för alla filer så man kan jämföra
dat <- read.table(paste("/home/Crisp/clinseq/BREASTv7/", files[2], sep=""), header = TRUE, as.is=TRUE)
# fread from data.table
plot(dat$segmented,pch=20)
dev.copy(jpeg,'breastv7_2.jpg')
dev.off()
head(dat)
dat10 <- subset(dat, dat$chromosome == "10")
plot(dat10$segmented,pch=20) 
head(dat10)

dat19cn <- subset(dat19, dat19$copynumber != "NA")  

savehistory(file="20160126")
loadhistory(file="20160126")
head(dat19cn)
bandplot(dat19$segmented)
  

  
  
  
#   list.files(path = ".", pattern = NULL, all.files = FALSE,
#              full.names = FALSE, recursive = FALSE,
#              ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
list.files(path="/home/Crisp/clinseq/BREASTv7/")
  