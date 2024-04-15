if (!require("genomation", quietly = TRUE)){
  install.packages("genomation")
}

# meth_diff_cut <- readRDS("data/meth_diff_cut_10_tiles_100.rds")

# Genome locations
# First load the annotation data; i.e the coordinates of promoters, TSS, intron and exons
ensembl_gtf2bed12_file <- "./annotations/Xtro10_ensemble_gtf2bed_ucsc_chr.bed"

if (!file.exists(ensembl_gtf2bed12_file)) {
  source("code/chromAlias_gtf2bed12.R")
}

ensembl_anot <- readTranscriptFeatures(ensembl_gtf2bed12_file, up.flank = 2000, down.flank = 2000) # 2000bp flank regions correspond to promoters. Upstream for + strand, downstream for - strand 

# Annotate hypermethylated CpGs ("target") with promoter/exon/intron
# information ("feature"). This function operates on GRanges objects, so we # first coerce the methylKit object to GRanges.
meth_diff_cut_annotated <- annotateWithGeneParts(target = as(meth_diff_cut,"GRanges"),
                                              feature = ensembl_anot, intersect.chr = TRUE, strand = TRUE)

gen_regions <- getTargetAnnotationStats(meth_diff_cut_annotated, percentage=FALSE,precedence=TRUE) %>% t %>% data.frame

# Summary of target set annotation
meth_diff_cut_annotated

meth_diff_cut_annotated_dir <- paste0("./data/meth_diff_cut_annotated",file_path_name,".rds")
saveRDS(meth_diff_cut_annotated, ascii=FALSE, file = meth_diff_cut_annotated_dir)

# View the distance to the nearest Transcription Start Site; the target.row column in the output indicates the row number in the initial target set
dist_tss <- getAssociationWithTSS(meth_diff_cut_annotated)
# hist(dist_tss$dist.to.feature)

# See whether the differentially methylated CpGs are within promoters,introns or exons; the order is the same as the target set
feature.sum <- genomation::getMembers(meth_diff_cut_annotated)
apply(feature.sum, 2, FUN = sum)

feature.dmr <- feature.sum %>% data.frame %>%
  mutate(location = case_when(
    prom == 1 ~ "promoter",
    exon == 1 ~ "exon",
    intron == 1 ~ "intron",
    TRUE ~ "intergenic"
  ))

# This can also be summarized for all differentially methylated CpGs
plot_genomic_location <- function() {genomation::plotTargetAnnotation(meth_diff_cut_annotated, main = "Genomic Features")}

tiff(paste0("./figures/genomic_location_", file_path_name, ".tiff"),
     width = 2000,
     height = 2000,
     units = "px",
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
)

plot_genomic_location()

dev.off()
