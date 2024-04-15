if (!require("circlize", quietly = TRUE)){
  install.packages("circlize")
}

if (!require("svglite", quietly = TRUE)){
  install.packages("svglite")
}

# meth_diff <- readRDS("./data/meth_diff_cut_10_tiles_100.rds")
# qvalue_cut <- 0.05
# meth_cut <- 10
# meth_diff_cut_hypo <- getMethylDiff(meth_diff,
#                                     difference = meth_cut,
#                                     qvalue = qvalue_cut,
#                                     type = "hypo")
# 
# meth_diff_cut_hyper <- getMethylDiff(meth_diff,
#                                      difference = meth_cut,
#                                      qvalue = qvalue_cut,
#                                      type = "hyper")
# feature_sum_annot <- readRDS("./data/annotated_df__cut_10_tiles100.rds")
# file_path_name <- "cut_10_tiles_100"


# prepare data frames with hypo and hyper methylated regions

DMR_hyper <- data.frame("chr" = meth_diff_cut_hyper$chr,
                        "start" = meth_diff_cut_hyper$start,
                        "end" = meth_diff_cut_hyper$end,
                        "meth.diff" = meth_diff_cut_hyper$meth.diff)

DMR_hypo <- data.frame("chr" = meth_diff_cut_hypo$chr,
                       "start" = meth_diff_cut_hypo$start,
                       "end" = meth_diff_cut_hypo$end,
                       "meth.diff" = meth_diff_cut_hypo$meth.diff)

bed_list <- list(DMR_hyper, DMR_hypo)

####All DMRs
# prepare data frame with all DMR with qvalue lower than the cut

all_DMRs <- data.frame(meth_diff) %>% filter(qvalue <= qvalue_cut)

DMR_all <- data.frame("chr" = all_DMRs$chr,
                      "start" = all_DMRs$start,
                      "end" = all_DMRs$end,
                      "meth.diff" = all_DMRs$meth.diff) 

DMR_all <- DMR_all %>%
  mutate(color = case_when(
    meth.diff >= 10 ~ "#E41A1C",
    meth.diff <= -10 ~ "#377EB8",
    meth.diff < 10 & meth.diff > -10 ~ "#D3D3D3"
  )) 




one <- c("clstn2","cacna2d3", "cacna1d", "cat2", "casr", "cacng3",
         "thbs4", "nox5", "mctp1", "eef2k", "melk", "cadps2", "ano1") #Calcium signalling
two <- c("minpp1", "aldh7a1", "tpi1" ,"eno3", "gckr", "uggt1") #Glucose metabolism
three <- c("tcf7l2", "ADCY5")
#four <- c("vti1a") #vesicle transport
#five <- c("igf1r") #pancreas development
four <- c("pnliprp2") #Lipase production
  

rlv.genes <- c(one, two, three, four)

# prepare data frame with selected genes that you want to highlight in the plot

hyper_color_rlv <- "#660000"
hypo_color_rlv <- "#000166"

annot <- feature_sum_annot %>% 
  filter(str_detect(external_gene_name, paste(rlv.genes, collapse = "|")))  %>% # Filter the names that you want to show in the plot
  tidyr::unite("loc", c("seqnames", "start", "end"), remove = FALSE) %>%
  dplyr::select(seqnames, start, end, mcols.meth.diff, external_gene_name, loc) %>%
  mutate(pointcolor = case_when(
    mcols.meth.diff >= 10 ~ hyper_color_rlv,
    mcols.meth.diff <= -10 ~ hypo_color_rlv
  )) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(one, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(one, collapse = "|"), 
                                                     "¹\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(two, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(two, collapse = "|"), 
                                                     "²\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(three, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(three, collapse = "|"), 
                                                     "³\\0"),
                                     external_gene_name)) %>%
  mutate(external_gene_name = ifelse(str_detect(external_gene_name, paste(four, collapse = "|")),
                                     str_replace_all(external_gene_name, 
                                                     paste(four, collapse = "|"), 
                                                     "⁴\\0"),
                                     external_gene_name))
  

annot$seqnames <- gsub("C", "c", annot$seqnames)

pointcolor <- annot %>% dplyr::select(pointcolor, loc)

DMR_all2 <- DMR_all
DMR_all2$loc <- paste(DMR_all$chr, DMR_all$start, DMR_all$end, sep = "_")

DMR_all2 <- left_join(DMR_all2, pointcolor, by = "loc") %>%
  mutate(color = if_else(is.na(pointcolor), color, pointcolor),
         cex = case_when(
           color == hyper_color_rlv ~ 1,
           color == hypo_color_rlv ~ 1,
           TRUE ~ 0.5
         )) %>%
  dplyr::select(-loc, -pointcolor)

DMR_all2 <- DMR_all2[order(DMR_all2$cex),]

circular_plot <- function() {
  circos.clear()
  
  circos.par(gap.after = c(rep(2, 9), 7), start.degree = 90)
  
  circos.initializeWithIdeogram(plotType = c("labels", "axis"),
                                species = "xenTro10"
  )
  circos.par("track.height" = 0.3)
  
  
  circos.genomicTrack(DMR_all2[,1:4], 
                      stack = FALSE, 
                      ylim = c(-100, 100),
                      panel.fun = function(region, value, ...) {
                        cex = DMR_all2[rownames(value), "cex"]
                        i = DMR_all2[rownames(value), "color"]
                        for(h in seq(-100, 100, by = 25)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
                        }
                        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, bg = i, ...)
                        circos.lines(CELL_META$cell.xlim, c(-10, -10), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(10, 10), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(0, 0), col = "#888888", lwd = 1, lty = "dotted")
                        circos.yaxis(side = "left", 
                                     at = c(-100, -50, 0, 50, 100),
                                     labels = paste0(c(-100, -50, 0, 50, 100), "%"),
                                     sector.index = get.all.sector.index()[1], 
                                     labels.cex = 0.3
                        )
                      })  
  circos.genomicLabels(annot, labels.column = 5, cex = 1.2, col= annot$pointcolor, line_lwd = 1, line_col="grey10",
                       side="inside", connection_height=0.025, labels_height=0.05, padding = 0.8)

}



tiff(paste0("./figures/circular_plot_", file_path_name, ".tiff"),
     width = 2000,
     height = 2000,
     units = "px",
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
)

circular_plot()

dev.off()

circular_plot()

svglite(paste0("./figures/circular_plot_", file_path_name, ".svg"),
        width = 16,
        height = 16,
        scaling = 2)

circular_plot()

dev.off()
