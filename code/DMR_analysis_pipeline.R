analysis_name <- "Differentially Methylated Regions"

source("code/read_bismark_coverage.R", echo = TRUE)

analysis_name <- "Differentially Methylated Regions"

analysis_object <- "CpG tiles"

myobj <- read_bismark_coverage(control_path, control_sample, lin_path, lin_sample)

params <- list(min_n = 4L, meth_cut = 10L, qvalue_cut = 0.05, 
               analysis = "tiles", tiling_window = "100", tiling_step = "100")

file_path_name <- paste0("cut_", params$meth_cut, "_", params$analysis,"_", params$tiling_window)

source("code/qc.R", echo = TRUE)

source("code/tiling_windows.R", echo = TRUE)

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
