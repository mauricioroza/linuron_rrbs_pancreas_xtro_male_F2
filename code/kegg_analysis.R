if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# Convert Ensembl gene to Entrezgene

res2 <- res %>%
  mutate(ensembl_gene_id = feature)

convert_ensembl_to_ens_df <- function(ensembl_genes_df) {
  
  # Convert ENSEMBL ID to Entrez ID using clusterProfiler::bitr
  
  ens_to_entrez <- bitr(ensembl_genes_df$feature, "ENSEMBL", "ENTREZID", org.xenTro.eg.db, drop = TRUE) %>%
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

KEGG <- function(genes, meth_status) {
  message(paste0("Analysing ", meth_status, " genes..."))
  
  # Convert Ensembl gene to Entrezgene
  
  entrez_df <- convert_ensembl_to_ens_df(genes)
  
  entrez <- entrez_df$ENTREZID %>% unique
  
  # KEGG gene enrichment analysis
  kegg_overrep <- enrichKEGG(gene         = entrez,
                         organism     = 'xtr',
                         pvalueCutoff = 0.5,
                         qvalueCutoff = 0.05,
                         pAdjustMethod = "BH",
                         use_internal_data = FALSE)
  
  if (is.null(kegg_overrep)) {
    print(paste0("No ", meth_status, " gene could be mapped"))
    assign(paste0("kegg_overrep_", meth_status), kegg_overrep, envir = .GlobalEnv)
    return(NULL)
  }
  
  
  
  # Convert Entrez gene number to gene names
  kegg_overrep <- setReadable(kegg_overrep, org.xenTro.eg.db, keyType = "ENTREZID")
  
  df_kegg_overrep <- data.frame(kegg_overrep) 
  
  if (nrow(df_kegg_overrep) == 0) {
    print(paste0("No significant KEGG pathways found for ", meth_status, " genes under qvalue cut = ", kegg_overrep@qvalueCutoff))
    return(NULL)
  } else {
    print(paste0(nrow(df_kegg_overrep), " significant KEGG pathways found for ", meth_status, " genes under qvalue cut = ", kegg_overrep@qvalueCutoff))
    
  }
  
  df_kegg_overrep_table <- df_kegg_overrep %>% datatable
  
  kegg_overrep_barplot <- kegg_overrep %>% 
    filter(., Count > 0) %>% 
    arrange(desc(Count)) %>% 
    barplot(., showCategory = 20, title = paste0("KEGG Pathways Enrichment ", analysis_name, " ", meth_status))

  # Save objects to local environment
  assign(paste0("kegg_overrep_", meth_status), kegg_overrep, envir = .GlobalEnv)
  assign(paste0("df_kegg_overrep_", meth_status), df_kegg_overrep, envir = .GlobalEnv)
  assign(paste0("df_kegg_overrep_table_", meth_status), df_kegg_overrep_table, envir = .GlobalEnv)
  assign(paste0("kegg_overrep_barplot_", meth_status), kegg_overrep_barplot, envir = .GlobalEnv)
  }

KEGG(res2, "all")

# Hypomethylated

hypo_genes <- res2 %>% dplyr::filter(mcols.meth.diff < -meth_cut)

KEGG(hypo_genes, "hypo")

# Hypermethylated

hyper_genes <- res2 %>% dplyr::filter(mcols.meth.diff > meth_cut)

KEGG(hyper_genes, "hyper")

# Promoters

feature_sum_annot2 <- feature_sum_annot %>% 
  mutate(ensembl_gene_id = feature)

prom_genes <- feature_sum_annot2 %>% dplyr::filter(prom == 1)
KEGG(prom_genes, "prom_all")

prom_hypo <- prom_genes %>% dplyr::filter(mcols.meth.diff < -meth_cut)
KEGG(prom_hypo, "prom_hypo")

prom_hyper <- prom_genes %>% dplyr::filter(mcols.meth.diff > meth_cut)
KEGG(prom_hyper, "prom_hyper")

#save RDS

save_RDS <- function(go_name) {
  if (exists(paste0(deparse(substitute(go_name)))) & !is.character(go_name)) {
    
    saveRDS(go_name, file = paste0("./data/", deparse(substitute(go_name)),"_", file_path_name, ".rds"))
    
  }
}

save_RDS(kegg_overrep_hyper)
save_RDS(kegg_overrep_hypo)


save_ggplot <- function (plot_name) {
  
  if (exists(paste0(deparse(substitute(plot_name))))) {
    
  ggsave(filename = paste0(deparse(substitute(plot_name)),"_", file_path_name, ".tiff"),
         path = "./figures/",
         plot = plot_name,
         device = "tiff",
         scale = 2, 
         dpi = 300,
         width = 15,
         height = 10,
         units = "cm",
         bg = "white")
  }
}

save_ggplot(kegg_overrep_barplot_all)
save_ggplot(kegg_overrep_barplot_hyper)
save_ggplot(kegg_overrep_barplot_hypo)

# save tables

save_table <- function(table_name) {
  
  if (exists(paste0(deparse(substitute(table_name))))) {
  
  write_xlsx(data.frame(table_name), path = paste0("./tables/", deparse(substitute(table_name)),"_", file_path_name, ".xlsx"))

  }
    
    }

save_table(df_kegg_overrep_all)
save_table(df_kegg_overrep_hyper)
save_table(df_kegg_overrep_hypo)
