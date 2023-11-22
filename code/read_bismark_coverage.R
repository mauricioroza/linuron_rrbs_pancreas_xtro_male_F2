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

# myobj <- read_bismark_coverage(control_path, control_sample, imz_path, imz_sample)
# myobj <- read_bismark_coverage(male_path, male_sample, female_path, female_sample)
# myobj <- read_bismark_coverage(control_male_path, control_male_sample, imz_male_path, imz_male_sample)


# Summarize methylated/unmethylated base counts over 100bp tilling windows accross genome
# myobj.tile <- tileMethylCounts(myobj.brain, win.size = 100, step.size = 100)