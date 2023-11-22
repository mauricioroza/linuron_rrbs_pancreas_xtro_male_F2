if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("methylKit", quietly = TRUE)){
  BiocManager::install("methylKit")
}

# Define the list containing the bismark coverage files.
file_path <- list.files("./methylation_coverage/", full.names = TRUE)
file_name <- list.files("./methylation_coverage/")
sample_name <- file_name %>% str_extract("\\d+_[A-Za-z]_[A-Za-z]\\d+")

# Sorting files
control_path <- file_path[grep("P_C", file_path)]
control_sample <- sample_name[grep("P_C", sample_name)]

lin_path <- file_path[grep("P_H", file_path)]
lin_sample <- sample_name[grep("P_H", sample_name)]


# Read the listed files into a methylRawList object making sure the other parameters are filled in correctly.
read_bismark_coverage <- function(group1_path, group1_names, group2_path, group2_names) {
  methRead(as.list(c(group1_path, group2_path)),
                    sample.id = as.list(c(group1_names, group2_names)),
                    pipeline = "bismarkCoverage",
                    assembly = "xenTro10",
                    treatment = c(rep(0, length(group1_names)), rep(1, length(group2_names))),
                    header = FALSE,
                    mincov = 3
  )
}

myobj <- read_bismark_coverage(control_path, control_sample, lin_path, lin_sample)

chr.alias <- read.table("annotations/Xtro_chromAlias.txt", header = TRUE)

chr.name1 <- as.list(chr.alias$genbank)
chr.name2 <- as.list(chr.alias$ucsc)

for(i in 1:length(myobj)){
  for(x in 1:length(chr.name1)){
    myobj[[i]]$chr[myobj[[i]]$chr == chr.name1[[x]]] <- chr.name2[[x]]
  }}
