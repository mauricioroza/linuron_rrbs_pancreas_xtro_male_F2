if (!require("circlize", quietly = TRUE)){
  install.packages("circlize")
}


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


# prepare data frame with selected genes that you want to highlight in the plot

annot <- feature_sum_annot %>% 
  filter(str_detect(description, "methyl"))  %>% # Filter the names that you want to show in the plot
  dplyr::select(seqnames, start, end, mcols.meth.diff, external_gene_name) %>%
  mutate(color = case_when(
    mcols.meth.diff >= meth_cut ~ "#E41A1C",
    mcols.meth.diff <= -meth_cut ~ "#377EB8"
  )
  )

annot$seqnames <- gsub("C", "c", annot$seqnames)

circular_plot <- function() {
  circos.clear()
  
  circos.par(gap.after = c(rep(2, 9), 7), start.degree = 90)
  
  circos.initializeWithIdeogram(plotType = c("labels", "axis"),
                                species = "xenTro10"
  )
  circos.par("track.height" = 0.3)
  
  circos.genomicTrack(DMR_all, 
                      stack = FALSE, 
                      ylim = c(-100, 100),
                      panel.fun = function(region, value, ...) {
                        cex = 0.5
                        i = ifelse(value[[1]] > meth_cut, "#E41A1C", ifelse(value[[1]] > -meth_cut, "#D3D3D3", "#377EB8"))
                        for(h in seq(-100, 100, by = 25)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
                        }
                        circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
                        circos.lines(CELL_META$cell.xlim, c(-meth_cut, -meth_cut), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(meth_cut, meth_cut), col = "#7D7D7D", lwd = 1, lty = "dashed")
                        circos.lines(CELL_META$cell.xlim, c(0, 0), col = "#888888", lwd = 1, lty = "dotted")
                        circos.yaxis(side = "left", 
                                     at = c(-100, -50, 0, 50, 100),
                                     labels = paste0(c(-100, -50, 0, 50, 100), "%"),
                                     sector.index = get.all.sector.index()[1], 
                                     labels.cex = 0.3
                        )
                      })  
  
  circos.par("track.height" = 0.15)
  
  circos.genomicDensity(bed_list, col = c("#FF000080", "#0000FF80"))
  
  # circos.genomicLabels(annot, labels.column=5, cex=0.5, col= annot$color, line_lwd=0.5, line_col="grey80", 
  #                      side="inside", connection_height=0.05, labels_height=0.04)
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
