analysis_name <- "Differentially Methylated Promoters"

source("code/read_bismark_coverage.R")
source("code/qc.R")
source("code/promoters_analysis.R")
source("code/filtering.R")
source("code/data_structure.R")
source("code/methylation_qvalue_cut_values.R")
source("code/differential_methylation.R")
source("code/genomic_location_DMR.R")
source("code/gene_annotation_overlap_DMR.R")
source("code/volcano_plot.R")
source("code/gene_ontology_analysis.R")
source("code/circular_plot.R")