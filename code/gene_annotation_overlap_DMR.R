if (!require("ChIPpeakAnno", quietly = TRUE)) {
  BiocManager::install("ChIPpeakAnno")
}

# Load mart and xenTro biomart table 
martfile <- c("./annotations/biomart_xenTro10.RData")
xen.biomart.file <- c("./annotations/xen.biomart.txt")

if (!file.exists(martfile)) {
  source("code/xenTro_biomart_table.R")
}

if (!file.exists(xen.biomart.file)) {
  source("code/xenTro_biomart_table.R")
}

xen.biomart <- read.table(xen.biomart.file, header = TRUE)
load(martfile)

# Table overlapping CpG regions with associated genes

# Add methylation difference and qvalues to the table
metadata <- getData(meth_diff_cut)

# Transform the differentially methylated regions into a GRanges object
meth_pos <- GRanges(
  seqnames = metadata$chr, 
  ranges = IRanges(start = metadata$start, end = metadata$end),
  mcols = cbind(metadata[, c(5:7)])
)

# Get annotation from Xenopus mart
xtro <- getAnnotation(mart, featureType = "TSS")

# Annotate DMRs with overlapping genes and promoters
res <- annotatePeakInBatch(meth_pos, AnnotationData = xtro, FeatureLocForDistance="TSS", output = "overlapping", maxgap = 2000, multiple = TRUE) %>% data.frame

# Remove scaffolds
res <- res %>% dplyr::filter(seqnames %in% paste0("chr", 1:10))

# Add gene names and description to the table
dmr_table_genes <- res %>%
  left_join(
    dplyr::select(xen.biomart, ensembl_gene_id, external_gene_name, description),
    by = c("feature" = "ensembl_gene_id"))

#dmr_table_genes[,c(1:3,6:8,10,14,18:19)] %>% datatable

# Annotate DMRs location (intron, exon or promoter)
anot <- res
anot$strand <- anot$feature_strand %>% replace_na(., "*")

anot <- annotateWithGeneParts(target = as(anot,"GRanges"),
                                    feature = ensembl_anot, intersect.chr = TRUE, strand = TRUE)
sum_anot <- data.frame(anot@members)
apply(sum_anot, 2, FUN = sum)

feature_sum_annot <- cbind(dmr_table_genes, sum_anot) %>%
  mutate(location = case_when(
    prom == 1 ~ "promoter",
    exon == 1 ~ "exon",
    intron == 1 ~ "intron",
    TRUE ~ "intergenic"
  ))
