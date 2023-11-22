analysis_name <- "Differentially Methylated Cytosines"

source("code/read_bismark_coverage.R")

analysis_name <- "Differentially Methylated Regions Male"

analysis_object <- "CpG sites"

myobj <- read_bismark_coverage(control_male_path, control_male_sample, imz_male_path, imz_male_sample)

params <- list(min_n = 4L, meth_cut = 10L, qvalue_cut = 0.05, 
               group = "male", analysis = "CpGs", tiling_window = "", tiling_step = "")

file_path_name <- paste0(params$group, "_cut_", params$meth_cut, "_", params$analysis, params$tiling_window)

source("code/qc.R")

min_n <- 4L

par(mfrow= c(1,1))

source("code/filtering.R")

source("code/data_structure.R")

qvalue_cut <- 0.05
meth_cut <- 10

source("code/differential_methylation.R")
source("code/genomic_location_DMR.R")
source("code/gene_annotation_overlap_DMR.R")
source("code/volcano_plot.R")
source("code/gene_ontology_analysis.R")
source("code/kegg_analysis.R")
