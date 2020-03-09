library(dplyr)

#Parse files to create a data frame with counts

files <- list.files(path = "significant_events/", pattern = "*.txt")
as_types <- c("a3ss", "a5ss", "mxe", "ri", "se")

files_aux <- gsub(pattern = ".txt", replacement = "", x = files)
files_aux <- gsub(pattern = "a3ss$|a5ss$|mxe$|ri$|se$", replacement = "", files_aux)

x <- table(files_aux)
tissues <- names(x)


counts <- rep(NA, length(files))
for (i in 1:length(files)) {
  #events <- read.table(paste0("./significant_events/", files[i]), sep = "\t", skip = 1)
  events <- readLines(paste0("./significant_events/", files[i]))
  counts[i] <- length(events) - 1 #header
}

ASE <- rep("NA", length(files))
for (i in 1:length(files)) {
  
  if(grepl("a3ss.txt$", files[i])){
    ASE[i] <- "A3SS"
  }
  if(grepl("a5ss.txt$", files[i])){
    ASE[i] <- "A5SS"
  }
  if(grepl("mxe.txt$", files[i])){
    ASE[i] <- "MXE"
  }
  if(grepl("se.txt$", files[i])){
    ASE[i] <- "SE"
  }
  if(grepl("ri.txt$", files[i])){
    ASE[i] <- "RI"
  }
  
}

data <- data.frame(Tissue = files_aux, ASE = ASE, Counts = counts)

head(data)

res <- data %>% group_by(Tissue) %>% 
    summarise(Total = sum(Counts)) %>%
    arrange(desc(Total)) %>%
    as.data.frame()

res

write.table(res, file = "Totals_by_tissue.tsv", sep = "\t", row.names = F)
write.table(data, file = "Significant_events.tsv", sep = "\t", row.names = F, quote = F)


