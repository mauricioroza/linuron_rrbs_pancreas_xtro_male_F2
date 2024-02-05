if (!require("genomation", quietly = TRUE)) {
  BiocManager::install("genomation")
}

if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# First load the annotation data; i.e the coordinates of promoters, TSS, intron and exons
ensembl_gtf2bed12_file <- "./annotations/Xtro10_ensemble_gtf2bed_ucsc_chr.bed"

if (!file.exists(ensembl_gtf2bed12_file)) {
  source("code/chromAlias_gtf2bed12.R")
}

martfile <- c("./annotations/biomart_xenTro10.RData")

if (!file.exists(martfile)) {
  source("code/xenTro_biomart_table.R")
}

load(martfile)

# ensembl_anot <- readTranscriptFeatures(ensembl_gtf2bed12_file, up.flank = 2000, down.flank = 0) # 2000bp flank regions correspond to promoters.
# 
# # Analysis of promoter regions
# 
# myobj <- regionCounts(myobj, c(ensembl_anot$promoters))
# 

##############

xentro_gene_list <- getBM(attributes = c("ensembl_gene_id", 'description', "chromosome_name", "start_position", "end_position", "strand"),
                      uniqueRows = TRUE,
                      mart = mart)

promoter_size <- 2000

promoters_df <- xentro_gene_list %>%
  dplyr::filter(chromosome_name %in% paste0(1:10)) %>%
  mutate(chromosome_name = paste0("chr", chromosome_name)) %>%
  mutate(strand = case_when(
    strand == "1" ~ "+",
    strand == "-1" ~ "-"
    )) %>%
  mutate(promoter_start = case_when(
    strand == "+" ~ start_position-promoter_size,
    strand == "-" ~ end_position
    )) %>%
  mutate(promoter_end = case_when(
    strand == "+" ~ start_position,
    strand == "-" ~ end_position+promoter_size
  ))

genes <- GRanges(
  seqnames = promoters_df$chromosome_name, 
  ranges = IRanges(start = promoters_df$start_position, end = promoters_df$end_position),
  strand = promoters_df$strand,
  mcols = cbind(promoters_df)
)

promoters <- GRanges(
  seqnames = promoters_df$chromosome_name, 
  ranges = IRanges(start = promoters_df$promoter_start, end = promoters_df$promoter_end),
  strand = promoters_df$strand,
  mcols = cbind(promoters_df)
)

myobj <- regionCounts(myobj, unique(promoters), cov.bases = 1, strand.aware = FALSE)
# myobj <- regionCounts(myobj, c(genes))
