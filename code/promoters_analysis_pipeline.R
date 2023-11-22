martfile <- c("./annotations/biomart_xenTro10.RData")

if (!file.exists(martfile)) {
  source("code/xenTro_biomart_table.R")
}

load(martfile)


analysis_name <- "Differentially Methylated Promoters"

source("code/read_bismark_coverage.R", echo = TRUE)

analysis_name <- "Differentially Methylated Promoters"

analysis_object <- "Promoters"


params <- list(min_n = 4L, meth_cut = 10L, qvalue_cut = 0.05, 
               analysis = "promoters", tiling_window = "", tiling_step = "")

file_path_name <- paste0("cut_", params$meth_cut, "_", params$analysis,"_", params$tiling_window)

source("code/qc.R", echo = TRUE)

source("code/promoters_analysis.R", echo = TRUE)

min_n <- params$min_n

source("code/filtering.R", echo = TRUE)

source("code/data_structure.R", echo = TRUE)

source("code/methylation_qvalue_cut_values.R", echo = TRUE)

source("code/differential_methylation.R", echo = TRUE)

source("code/genomic_location_DMR.R", echo = TRUE)

source("code/gene_annotation_overlap_DMR.R", echo = TRUE)

source("code/volcano_plot.R", echo = TRUE)

source("code/gene_ontology_analysis.R", echo = TRUE)

source("code/kegg_analysis.R", echo = TRUE)

