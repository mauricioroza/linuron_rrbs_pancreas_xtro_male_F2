# Select dataframe for gene list for gene set enrichment analysis


annotate_dmr <- function (meth_diff_cut) {
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
  
  # Add gene names and description to the table
  dmr_table_genes <- res %>%
    left_join(
      dplyr::select(xen.biomart, ensembl_gene_id, external_gene_name, description),
      by = c("feature" = "ensembl_gene_id"))
  return(dmr_table_genes)
}

dmr_table_genes <- annotate_dmr(meth_diff)

# result <- res %>%
#   group_by(feature) %>%
#   summarise(mcols.meth.diff = mean(mcols.meth.diff, na.rm = TRUE))
# 
# unique_rows <- feature_sum_annot[!duplicated(feature_sum_annot$feature), ]
# rownames(unique_rows) <- NULL


gsea_list_df <- dmr_table_genes %>%
  dplyr::select(feature, mcols.meth.diff, mcols.pvalue, mcols.qvalue) %>%
  mutate(metric = mcols.meth.diff*(-log10(mcols.qvalue))) %>%
  mutate(ensembl_gene_id = feature)

unique_rows <- gsea_list_df[!duplicated(gsea_list_df$feature), ]
rownames(unique_rows) <- NULL

gsea_list_df <- unique_rows





#  column_to_rownames("ensembl_gene_id")

# Select metric for ranking the genes and arrange in decreasing order
gsea_list <- gsea_list_df$metric
names(gsea_list) <- gsea_list_df$feature
gsea_list <- gsea_list[order(gsea_list, decreasing = TRUE)]

# Gene set enrichment analysis for GO terms
gseaGO <- gseGO(gsea_list,
                ont = "BP",
                OrgDb = org.xenTro.eg.db,
                keyType = "ENSEMBL",
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH"
)


###########################################################################

convert_ensembl_to_ens_df <- function(ensembl_genes_df) {
  
  # Convert ENSEMBL ID to Entrez ID using clusterProfiler::bitr
  
  ens_to_entrez <- bitr(ensembl_genes_df$ensembl_gene_id, "ENSEMBL", "ENTREZID", org.xenTro.eg.db, drop = TRUE) %>%
    left_join(ensembl_genes_df, by = c("ENSEMBL" = "ensembl_gene_id"))
  
  # Get gene names for Entrez IDs and join in the df
  
  entrez_names <- getBM(attributes = c("entrezgene_accession", 'entrezgene_id', "entrezgene_description"),
                        filters = 'entrezgene_id',
                        values = ens_to_entrez$ENTREZID,
                        uniqueRows = TRUE,
                        mart = mart) %>%
    distinct()
  
  entrez_names$entrezgene_id <- as.character(entrez_names$entrezgene_id)
  
  ens_to_entrez <- ens_to_entrez %>%
    left_join(entrez_names, by = c("ENTREZID" = "entrezgene_id"))
  
  return(ens_to_entrez)
}


# Use function on gene list data frame
gsea_list_df_entrez <- convert_ensembl_to_ens_df(gsea_list_df)

# Select metric to rank the genes  
gsea_list_entrez <- gsea_list_df_entrez$metric
names(gsea_list_entrez) <- gsea_list_df_entrez$ENTREZID

# Arrange list in decreasing order
gsea_list_entrez <- gsea_list_entrez[order(gsea_list_entrez, decreasing = TRUE)]

# Run Gene Set Enrichment Analysis for KEGG Pathways
kegg_gsea <- gseKEGG(geneList     = gsea_list_entrez,
                     organism     = 'xtr',
                     minGSSize    = 10,
                     pvalueCutoff = 0.05,
                     verbose      = FALSE)



