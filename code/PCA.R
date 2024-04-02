library(tidyverse)
library(readxl)
library(ChIPpeakAnno)
library(methylKit)

######################################################################

perc_meth_table <- function(meth.diff, perc.meth) {
  
  mydiff <- getMethylDiff(meth.diff,
                          difference=meth.cut,
                          qvalue=qvalue.cut) %>% getData()
  
  perc.meth.dm <- perc.meth[rownames(perc.meth) %in% rownames(mydiff),]
  
  perc.meth.dm.merged <- cbind(mydiff, perc.meth.dm)
  
  meth_pos.pca <- GRanges(
    seqnames = perc.meth.dm.merged$chr, 
    ranges = IRanges(start = perc.meth.dm.merged$start, end = perc.meth.dm.merged$end),
    mcols = perc.meth.dm.merged[, c(5:ncol(perc.meth.dm.merged))]
  )
  
  martfile <- c("./data/mart_xtropicalis_gene_ensembl.RData")
  
  load(martfile)
  
  xtro <- getAnnotation(mart, featureType = "TSS")
  
  pca <- annotatePeakInBatch(meth_pos.pca, AnnotationData = xtro, FeatureLocForDistance="TSS", output = "overlapping", maxgap = 2000, multiple = TRUE) %>% data.frame
  
  xen.biomart <- biomaRt::getBM(attributes = 
                                  c("ensembl_gene_id","external_gene_name","chromosome_name", 
                                    "start_position","end_position","description"
                                  ),       
                                filters    = "",
                                values = "",
                                mart       = mart) 
  
  t.pca <- pca %>%
    left_join(
      dplyr::select(xen.biomart, ensembl_gene_id, external_gene_name, description),
      by = c("feature" = "ensembl_gene_id"))
  return(t.pca)
}

rlv_genes_table <- function(t.pca, rlv.genes, phenotype) {
  
  rlv.genes.table <- t.pca %>%
    filter(str_detect(external_gene_name, paste(rlv.genes, collapse = "|")))
  
  t.rlv.genes.table <- rlv.genes.table %>% t %>% data.frame
  
  
  
  pca.rlv.genes.table <- t.rlv.genes.table[grep("", rownames(t.rlv.genes.table)), ] %>% 
    rownames_to_column('ID') %>%
    mutate_at('ID', str_replace, "mcols.", "")
  
  colnames(pca.rlv.genes.table) <- c("ID", t.rlv.genes.table["external_gene_name",])
  
  rlv.table.pca <- merge(phenotype, pca.rlv.genes.table, by = "ID")
  
  return(rlv.table.pca)
}

PCA_phenotype <- function(rlv.table.pca, plot_title) {
  
  rlv.table.pca[,c(3:ncol(rlv.table.pca))] <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% mutate_if(is.character,as.numeric)
  
  library(missMDA)
  imp<-imputePCA(rlv.table.pca, 
                 
                 quali.sup = 1:2)
  imp.prop.variables<-imp$completeObs
  
  library(factoextra)
  library(FactoMineR)
  pca.prop <- PCA(imp.prop.variables, 
                  scale.unit = TRUE, 
                  graph = FALSE,
                  quali.sup = 1:2)
  
  
  ##
  biplot.genexpression<- fviz_pca_biplot(pca.prop, repel = TRUE,
                                         axes = c(1, 2),
                                         label = "var", habillage="treatment", addEllipses=FALSE, 
                                         col.var = "Gray31", # Variables color
                                         col.ind = c("Black"), 
                                         title = plot_title
  )
  
  
  assign(paste0("pca_", plot_title), pca.prop, envir = .GlobalEnv)
  biplot.genexpression

}


unite <- readRDS("./data/meth_united_cut_10_tiles_100.rds")

pm <- percMethylation(unite) %>% 
  data.frame
colnames(pm) <- str_replace(colnames(pm), "X\\d{1,2}_P_", "")

pancreas.meth.dir <- c("./data/meth_diff_cut_10_tiles_100.rds")
pancreas.meth <- readRDS(pancreas.meth.dir)

meth.cut <- 10
qvalue.cut <- 0.05

perc_meth_table <- perc_meth_table(pancreas.meth, pm)

fat <- read_excel("data/fat_body_fatty_acids.xlsx") %>%
  dplyr::select(!c('20:3n3', "20:4n6", "DPA", "DHA"))

liver <- read_excel("data/liver_fatty_acids.xlsx") %>%
  dplyr::select(!c('20:3n3', "18:3n3"))

phenotype <- read_excel("data/phenotype_data.xlsx")

# phenotype <- phenotype %>%
#   dplyr::select(treatment, ID, body_weight, cholesterol, triglycerids, glucose)

phenotype <- phenotype %>%
  dplyr::select(treatment, ID, glucose, body_weight, cholesterol, triglycerids)

one <- c("pnliprp2") #Lipase production
two <- c("minpp1", "aldh7a1", "tpi1" ,"eno3", "gckr", "uggt1") #Glucose metabolism
three <- c("clstn2","cacna2d3", "cacna1d", "cat2", "casr", "cacng3",
           "thbs4", "nox5", "mctp1", "eef2k", "melk", "cadps2", "ano1") #Calcium signalling
four <- c("vti1a") #vesicle transport
five <- c("igf1r") #pancreas development
six <- c("tcf7l2", "ADCY5", "gckr") #

phen <- merge(fat, phenotype, by = c("ID", "treatment"))

gene_list <- c(one, two, three, six)
rlv.table.pca <- rlv_genes_table(perc_meth_table, gene_list, phenotype)
PCA_phenotype(rlv.table.pca,  "pca_all_genes")

phenotype.names <- names(phen)
pca.var.names <- names(rlv.table.pca)[-c(1:2)]

df.groups <- data.frame(var = pca.var.names) 

df.groups <- df.groups %>%
  mutate(group = case_when(
    str_detect(var, paste(phenotype.names, collapse = "|")) ~ "Phenotype",
    str_detect(var, paste(one, collapse = "|")) ~ "Lipase production",
    str_detect(var, paste(two, collapse = "|")) ~ "Glucose metabolism",
    str_detect(var, paste(three, collapse = "|")) ~ "Calcium signalling",
    str_detect(var, paste(six, collapse = "|")) ~ "Insulin secretion"
  ))

gene_func <- as.factor(df.groups$group)
names(gene_func) <- df.groups$var

library(RColorBrewer)
col <- brewer.pal(n = length(group), "Set1")


var_plot <- fviz_pca_var(pca_pca_all_genes, repel = TRUE, col.var = gene_func, palette = col)

ind_plot <- fviz_pca_ind(pca_pca_all_genes, 
                         habillage=as.factor(rlv.table.pca$treatment),
                         label = "var")

biplot <- fviz_pca_biplot(pca_pca_all_genes,
                          fill.ind = as.factor(rlv.table.pca$treatment),
                          col.ind = "black",
                          pointshape = 21,
                          palette = col,
                          col.var = gene_func,
                          label = "var",
                          repel = TRUE) +
  geom_point(shape = as.factor(rlv.table.pca$treatment)) +
  scale_fill_manual(values = c("black", "red")) +
  labs(fill = "Treatment", color = "Gene Function")  # Change legend title
  
  
  
# Contributions of variables to PC1
fviz_contrib(pca_pca_all_genes, choice = "var", axes = 1)
# Contributions of variables to PC2
fviz_contrib(pca_pca_all_genes, choice = "var", axes = 2)

fviz_dend(pca_pca_all_genes)


feature_sum_annot <- readRDS("./data/annotated_df__cut_10_tiles100.rds")
rlv_genes_table_print <- feature_sum_annot %>%
  filter(external_gene_name %in% c(one,two,three)) %>%
  arrange(match(external_gene_name, c(one,two,three))) %>%
  dplyr::select(external_gene_name, location, mcols.meth.diff, mcols.qvalue, description) %>%
  mutate(mcols.meth.diff = round(mcols.meth.diff, digits = 1),
         mcols.qvalue = round(mcols.qvalue, digits = 3),
         description = sub("\\s*\\[Source.*", "", description)
  )
clipr::write_clip(rlv_genes_table_print, object_type = "table")
