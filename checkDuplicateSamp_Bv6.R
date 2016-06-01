



samples_no_T <- gsub("T.*panel_v1","",samples)

for (i in 1:length(samples_no_T)) {
  dupl <- which(samples_no_T[i]==samples_no_T)
  print(paste("Sample", i, samples_no_T[i], "exists", length(dupl), "times:", samples[dupl]))
}


