library(tidyverse)
library(rstatix)

qc <- read.csv("./data/multiqc_general_statistics.csv")

qc.p <- qc %>%
  filter(str_detect(Sample.Name, "_P_")) %>%
  mutate(across(starts_with("X.."), parse_number))

qc.table <- qc.p %>%
  get_summary_stats(type = "mean_sd")

qc.table

