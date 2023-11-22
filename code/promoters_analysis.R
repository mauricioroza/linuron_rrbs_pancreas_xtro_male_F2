if (!require("genomation", quietly = TRUE)) {
  BiocManager::install("genomation")
}

# First load the annotation data; i.e the coordinates of promoters, TSS, intron and exons
ensembl_gtf2bed12_file <- "./annotations/Xtro10_ensemble_gtf2bed_assembly_chr.bed"

if (!file.exists(ensembl_gtf2bed12_file)) {
  source("code/chromAlias_gtf2bed12.R")
}

ensembl_anot <- readTranscriptFeatures(ensembl_gtf2bed12_file, up.flank = 2000, down.flank = 2000) # 2000bp flank regions correspond to promoters. Upstream for + strand, downstream for - strand 

# Analysis of promoter regions

myobj <- regionCounts(myobj, c(ensembl_anot$promoters))