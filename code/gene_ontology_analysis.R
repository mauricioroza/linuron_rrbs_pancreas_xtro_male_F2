if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!require("org.xenTro.eg.db", quietly = TRUE)) {
  source("code/org.db_package.R")
}

if (!require("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

if (!require("DT", quietly = TRUE)){
  install.packages("DT")
}

gene_ontology <- function(gene_list, meth_status) {

  message(paste0("Analysing ", meth_status, " genes..."))
  
  # Remove duplicated genes
  gene_list <- gene_list %>% na.omit %>% unique
  
  #gene ontology terms BP biological function
  ggo <- groupGO(gene     = gene_list,
                       OrgDb    = org.xenTro.eg.db,
                       keyType = "ENSEMBL",
                       ont      = "BP",
                       readable = TRUE,
                       level = 2)
  ggo_result <- ggo@result
  
  ggo_table <- ggo@result %>% datatable
  
  ggo_plot <- ggo %>% filter(., Count > 0) %>% arrange(desc(Count)) %>% barplot(., showCategory = 20, title = paste0("GO BP ", analysis_name, " ", meth_status))
  
  #gene onlogy terms CC cellular compartment
  ggo_cc <- groupGO(gene     = gene_list,
                          OrgDb    = org.xenTro.eg.db,
                          keyType = "ENSEMBL",
                          ont      = "CC",
                          readable = TRUE)
  
  ggo_cc_table <- ggo_cc@result %>% datatable
  
  ggo_cc_barplot <- ggo_cc %>% filter(., Count > 0) %>% arrange(desc(Count)) %>% barplot(., showCategory = 20)
  
  #gene onlogy terms MF molecular function
  ggo_mf <- groupGO(gene     = gene_list,
                          OrgDb    = org.xenTro.eg.db,
                          keyType = "ENSEMBL",
                          ont      = "MF",
                          readable = TRUE)
  
  ggo_mf_table <- ggo_mf@result %>% datatable
  
  ggo_mf_barplot <- ggo_mf %>% filter(., Count > 0) %>% arrange(desc(Count)) %>% barplot(., showCategory = 20)
  
  #GO gene set over-representation
  go_overrep <- enrichGO(gene          = gene_list,
                        OrgDb         = org.xenTro.eg.db,
                        keyType = "ENSEMBL",
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 0.05,
                        readable      = TRUE)
  
  go_overrep_df <- data.frame(go_overrep)
  
  if (nrow(go_overrep_df) == 0) {
    print(paste0("No significant GO terms over-represented for ", meth_status, " genes under qvalue cut = ", go_overrep@qvalueCutoff))
    
    go_overrep <- paste0("No significant GO terms over-represented")
    go_overrep_table <- paste0("No significant GO terms over-represented")
    go_overrep_barplot <- paste0("No significant GO terms over-represented")
    
    assign(paste0("go_overrep_", meth_status), go_overrep, envir = .GlobalEnv)
    assign(paste0("go_overrep_table_", meth_status), go_overrep_table, envir = .GlobalEnv)
    assign(paste0("go_overrep_barplot_", meth_status), go_overrep_barplot, envir = .GlobalEnv)
    
  } else {
    print(paste0(nrow(go_overrep_df), " significant GO terms over-represented found for ", meth_status, " genes under qvalue cut = ", go_overrep@qvalueCutoff))
    
    go_overrep_table <- data.frame(go_overrep) %>% datatable
    
    go_overrep_barplot <- barplot(go_overrep, color = "qvalue", showCategory = 10, title = paste0("GO Over-representation ", analysis_name, " ", meth_status))
    go_overrep_barplot
    
    
    assign(paste0("go_overrep_", meth_status), go_overrep, envir = .GlobalEnv)
    assign(paste0("go_overrep_table_", meth_status), go_overrep_table, envir = .GlobalEnv)
    assign(paste0("go_overrep_barplot_", meth_status), go_overrep_barplot, envir = .GlobalEnv)
  }
  
  print(meth_status)

  # Save objects to local environment
  assign(paste0("ggo_plot_", meth_status), ggo_plot, envir = .GlobalEnv)
  assign(paste0("ggo_result_", meth_status), ggo_result, envir = .GlobalEnv)
  assign(paste0("ggo_table_", meth_status), ggo_table, envir = .GlobalEnv)

  }

# All

gene_list <- res

gene_ontology(gene_list$feature, "all")

# Hypo

hypo_genes <- res %>% dplyr::filter(mcols.meth.diff < -meth_cut)

gene_ontology(hypo_genes$feature, "hypo")

# Hyper

hyper_genes <- res %>% dplyr::filter(mcols.meth.diff > meth_cut)

gene_ontology(hyper_genes$feature, "hyper")

# Promoters

prom_genes <- feature_sum_annot %>% dplyr::filter(prom == 1)
gene_ontology(prom_genes$feature, "prom_all")

prom_hypo <- prom_genes %>% dplyr::filter(mcols.meth.diff < -meth_cut)
gene_ontology(prom_hypo$feature, "prom_hypo")

prom_hyper <- prom_genes %>% dplyr::filter(mcols.meth.diff > meth_cut)
gene_ontology(prom_hyper$feature, "prom_hyper")

# save plots 

save_ggplot <- function (plot_name) {
  
  if (exists(paste0(deparse(substitute(plot_name)))) & !is.character(plot_name)) {
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

save_ggplot(ggo_plot_all)
save_ggplot(go_overrep_barplot_all)
save_ggplot(go_overrep_barplot_hyper)
save_ggplot(go_overrep_barplot_hypo)
save_ggplot(go_overrep_barplot_prom_hypo)


save_RDS <- function(go_name) {
  if (exists(paste0(deparse(substitute(go_name)))) & !is.character(go_name)) {
    
    saveRDS(go_name, file = paste0("./data/", deparse(substitute(go_name)),"_", file_path_name, ".rds"))
    
  }
}

save_RDS(ggo_result_all)
save_RDS(go_overrep_all)
save_RDS(go_overrep_hyper)
save_RDS(go_overrep_hypo)
save_RDS(go_overrep_prom_hypo)
save_RDS(go_overrep_prom_hyper)


# save tables

save_table <- function(table_name) {
  if (exists(paste0(deparse(substitute(table_name)))) & !is.character(table_name)) {
    
    write_xlsx(data.frame(table_name), path = paste0("./tables/", deparse(substitute(table_name)),"_", file_path_name, ".xlsx"))
    
  }
}

save_table(ggo_result_all)
save_table(go_overrep_all)
save_table(go_overrep_hyper)
save_table(go_overrep_hypo)
save_table(go_overrep_prom_hypo)
save_table(go_overrep_prom_hyper)

