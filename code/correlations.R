library(tidyverse)
library(readxl)
library(methylKit)
library(ChIPpeakAnno)

#' correlation_matrix
#' Creates a publication-ready / formatted correlation matrix, using `Hmisc::rcorr` in the backend.
#'
#' @param df dataframe; containing numeric and/or logical columns to calculate correlations for
#' @param type character; specifies the type of correlations to compute; gets passed to `Hmisc::rcorr`; options are `"pearson"` or `"spearman"`; defaults to `"pearson"`
#' @param digits integer/double; number of decimals to show in the correlation matrix; gets passed to `formatC`; defaults to `3`
#' @param decimal.mark character; which decimal.mark to use; gets passed to `formatC`; defaults to `.`
#' @param use character; which part of the correlation matrix to display; options are `"all"`, `"upper"`, `"lower"`; defaults to `"all"`
#' @param show_significance boolean; whether to add `*` to represent the significance levels for the correlations; defaults to `TRUE`
#' @param replace_diagonal boolean; whether to replace the correlations on the diagonal; defaults to `FALSE`
#' @param replacement character; what to replace the diagonal and/or upper/lower triangles with; defaults to `""` (empty string)
#'
#' @return a correlation matrix
#' @export
#'
#' @examples
#' `correlation_matrix(iris)`
#' `correlation_matrix(mtcars)`
correlation_matrix <- function(df, 
                               type = "spearman",
                               digits = 3, 
                               decimal.mark = ".",
                               use = "all", 
                               show_significance = TRUE, 
                               replace_diagonal = FALSE, 
                               replacement = ""){
  
  # check arguments
  stopifnot({
    is.numeric(digits)
    digits >= 0
    use %in% c("all", "upper", "lower")
    is.logical(replace_diagonal)
    is.logical(show_significance)
    is.character(replacement)
  })
  # we need the Hmisc package for this
  require(Hmisc)
  
  # retain only numeric and boolean columns
  isNumericOrBoolean = vapply(df, function(x) is.numeric(x) | is.logical(x), logical(1))
  if (sum(!isNumericOrBoolean) > 0) {
    cat('Dropping non-numeric/-boolean column(s):', paste(names(isNumericOrBoolean)[!isNumericOrBoolean], collapse = ', '), '\n\n')
  }
  df = df[isNumericOrBoolean]
  
  # transform input data frame to matrix
  x <- as.matrix(df)
  
  # run correlation analysis using Hmisc package
  correlation_matrix <- Hmisc::rcorr(x, type = type)
  R <- correlation_matrix$r # Matrix of correlation coeficients
  p <- correlation_matrix$P # Matrix of p-value 
  
  # transform correlations to specific character format
  Rformatted = formatC(R, format = 'f', digits = digits, decimal.mark = decimal.mark)
  
  # if there are any negative numbers, we want to put a space before the positives to align all
  if (sum(!is.na(R) & R < 0) > 0) {
    Rformatted = ifelse(R > 0, paste0(' ', Rformatted), Rformatted)
  }
  
  # add significance levels if desired
  if (show_significance) {
    # define notions for significance levels; spacing is important.
    stars <- ifelse(is.na(p), "   ", ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "*  ", "   "))))
    Rformatted = paste0(Rformatted, stars)
  }
  # build a new matrix that includes the formatted correlations and their significance stars
  Rnew <- matrix(Rformatted, ncol = ncol(x))
  rownames(Rnew) <- colnames(x)
  colnames(Rnew) <- paste(colnames(x), "", sep =" ")
  
  # replace undesired values
  if (use == 'upper') {
    Rnew[lower.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (use == 'lower') {
    Rnew[upper.tri(Rnew, diag = replace_diagonal)] <- replacement
  } else if (replace_diagonal) {
    diag(Rnew) <- replacement
  }
  
  return(Rnew)
}


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
  
  rlv.table.pca[,c(3:ncol(rlv.table.pca))] <- rlv.table.pca[,c(3:ncol(rlv.table.pca))] %>% 
    mutate_if(is.character,as.numeric)
  return(rlv.table.pca)
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

liver <- read_excel("data/liver_fatty_acids.xlsx") %>%
  dplyr::select(!c('20:3n3', "18:3n3"))

phenotype <- read_excel("data/phenotype_data.xlsx")

phenotype <- phenotype %>%
  dplyr::select(treatment, ID, body_weight, cholesterol, triglycerids, glucose)

one <- c("pnliprp2", "lmf1") #Lipase production
two <- c("minpp1", "aldh7a1", "tpi1" ,"eno3", "gckr") #Glucose metabolism
three <- c("clstn2","cacna2d3", "thbs4", "nox5", "mctp1", "eef2k", "melk", "cadps2", "ano1") #Calcium signalling
four <- c("vti1a") #vesicle transport
five <- c("igf1r") #pancreas development
six <- c("tcf7l2", "ADCY5")

phen <- merge(fat, phenotype, by = c("ID", "treatment"))

df <- rlv_genes_table(perc_meth_table, c(one, two, three, four, five), phen)

corr_table <- correlation_matrix(df, type = "spearman")

corr_table <- corr_table[-c((ncol(phen)-1):ncol(corr_table)), -c(1:(ncol(phen)-2))] %>% data.frame %>%
  rownames_to_column(var = "Gene")


phen.liver <- merge(liver, phenotype, by = c("ID", "treatment"))

df.liver <- rlv_genes_table(perc_meth_table, c(one, two, three, four, five,six), phen.liver)

corr_table.liver <- correlation_matrix(df.liver, type = "spearman")

corr_table.liver <- corr_table.liver[-c((ncol(phen.liver)-1):ncol(corr_table.liver)), -c(1:(ncol(phen.liver)-2))] %>% data.frame %>%
  rownames_to_column(var = "Gene")
