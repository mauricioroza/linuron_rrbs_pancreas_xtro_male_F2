if (!require("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}

# Save biomaRt table

# Select a mart and data set        
mart <- biomaRt::useDataset(dataset = "xtropicalis_gene_ensembl",         
                            mart    = useMart("ENSEMBL_MART_ENSEMBL",       
                                              host    = "https://www.ensembl.org"))   

# Save mart file
martfile <- c("./annotations/biomart_xenTro10.RData")
save(mart, ascii=FALSE, file=martfile)
load(martfile)
mart

#Download annotations from ensembl biomaRt
xen.biomart <- biomaRt::getBM(attributes = 
                                c("ensembl_gene_id","external_gene_name","chromosome_name", 
                                  "start_position","end_position","description"
                                ),       
                              filters    = "",
                              values = "",
                              mart       = mart) 


#Change chomosome names
chromAlias_path <- "./annotations/xenTro10chromAlias.txt.gz"

if (!file.exists("./annotations/xenTro10chromAlias.txt.gz")) {
  chromAlias_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/xenTro10/database/chromAlias.txt.gz"
  download.file(chromAlias_url, destfile = chromAlias_path)
}

chr_alias <- read.table(gzfile(chromAlias_path), header = FALSE, sep = "\t")

## add a column to chomAlias with Ensembl chr name pattern that is found in gtf file

chr_alias <- chr_alias %>% 
  pivot_wider(names_from = V3, values_from = V1) %>% 
  rename(ucsc = V2) %>%
  mutate(ensembl = genbank) %>%
  mutate(ensembl = ifelse(row_number() %in% 1:11, c(1, 10, 2, 3, 4, 5, 6, 7, 8, 9, "MT"), ensembl))

saveRDS(chr_alias, "./annotations/Xtro_chromAlias.txt")

xen.biomart = xen.biomart %>%
  ## join the relevant alias column into xen biomart
  left_join(
    dplyr::select(chr_alias, ensembl, assembly),
    by = c("chromosome_name" = "ensembl")
  ) %>%
  ## replace all chromosome_names with ucsc value (if not NA)
  mutate(chromosome_name = coalesce(assembly, chromosome_name)) %>%
  ## drop ucsc columns
  dplyr::select(-assembly)

write.table(xen.biomart, file = "./annotations/xen.biomart.txt", row.names = FALSE, col.names = TRUE)
