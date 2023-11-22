# Load gtf that was converted to bed12 file

Xtro10_ensemble_gtf2bed <- read.table("./annotations/Xtro10_ensemble_gtf2bed.bed", header = FALSE)
levels(as.factor(Xtro10_ensemble_gtf2bed$V1))

# Change chromosome names to ucsc format

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

write.table(chr_alias, file = "./annotations/Xtro_chromAlias.txt", 
            row.names = FALSE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)


Xtro10_ensemble_gtf2bed_ucsc_Chr = data.frame(Xtro10_ensemble_gtf2bed) %>%
  left_join(
    dplyr::select(chr_alias, ensembl, ucsc),
    by = c("V1" = "ensembl")
  ) %>%
  mutate(V1 = coalesce(ucsc, V1)) %>%
  dplyr::select(-ucsc)

levels(as.factor(Xtro10_ensemble_gtf2bed_ucsc_Chr$V1))

#Save .bed file with changed chromosome name format
write.table(Xtro10_ensemble_gtf2bed_ucsc_Chr, file = "./annotations/Xtro10_ensemble_gtf2bed_ucsc_chr.bed", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

Xtro10_ensemble_gtf2bed_ucsc_Chr <- read.table("./annotations/Xtro10_ensemble_gtf2bed_ucsc_chr.bed", header = FALSE)
