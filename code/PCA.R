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
  
  variables.pca <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% mutate_if(is.character,as.numeric)
  
  library(missMDA)
  imp<-imputePCA(variables.pca)
  imp.prop.variables<-imp$completeObs
  
  library(factoextra)
  library(FactoMineR)
  pca.prop <- PCA(imp.prop.variables, scale.unit = TRUE, graph = FALSE)
  
  
  ##
  biplot.genexpression<- fviz_pca_biplot(pca.prop, repel = TRUE,
                                         axes = c(1, 2),
                                         label = "var", habillage=as.factor(rlv.table.pca$treatment), addEllipses=FALSE, 
                                         col.var = "Gray31", # Variables color
                                         col.ind = c("Black"), 
                                         title = plot_title
  )
  
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
  dplyr::select(!c('20:3n3', "20:4n6", "DPA"))

phenotype <- read_excel("data/phenotype_data.xlsx")

phenotype <- phenotype %>%
  dplyr::select(treatment, ID, body_weight, cholesterol, triglycerids, glucose)

one <- c("pnliprp2", "lmf1") #Lipase production
two <- c("minpp1", "aldh7a1", "tpi1" ,"eno3", "gckr") #Glucose metabolism
three <- c("clstn2","cacna2d3", "thbs4", "nox5", "mctp1", "eef2k", "melk", "cadps2") #Calcium signalling
four <- c("vti1a") #vesicle transport
five <- c("igf1r") #pancreas development

phen <- merge(fat, phenotype, by = c("ID", "treatment"))

rlv.table.pca <- rlv_genes_table(perc_meth_table, c(one), phen)

PCA_phenotype(rlv.table.pca,  "PCA Lipases Genes Brain")
