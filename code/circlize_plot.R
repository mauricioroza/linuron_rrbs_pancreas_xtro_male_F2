library(circlize)

DMR_hyper <- data.frame("chr" = myDiff10p.hyper.pancreas$chr,
                        "start" = myDiff10p.hyper.pancreas$start,
                        "end" = myDiff10p.hyper.pancreas$end,
                        "meth.diff" = myDiff10p.hyper.pancreas$meth.diff)

DMR_hypo <- data.frame("chr" = myDiff10p.hypo.pancreas$chr,
                        "start" = myDiff10p.hypo.pancreas$start,
                        "end" = myDiff10p.hypo.pancreas$end,
                        "meth.diff" = myDiff10p.hypo.pancreas$meth.diff)

bed_list <- list(DMR_hyper, DMR_hypo)


####All DMRs

all.DMRs <- data.frame(myDiff.pancreas) %>% filter(qvalue <= 0.05)

DMR_all <- data.frame("chr" = all.DMRs$chr,
                       "start" = all.DMRs$start,
                       "end" = all.DMRs$end,
                       "meth.diff" = all.DMRs$meth.diff)

annot <- t.pancreas %>% 
  filter(str_detect(description, "methyl")) %>%
  dplyr::select(seqnames, start, end, mcols.meth.diff, external_gene_name) %>%
  mutate(color = case_when(
    mcols.meth.diff >= 10 ~ "#E41A1C",
    mcols.meth.diff <= -10 ~ "#377EB8"
  )
         )


tiff("./figures/circos.tiff", 
     width = 2000, 
     height = 2000, 
     units = "px", 
     pointsize = 12,
     compression = "none",
     bg = "white",
     res = 300
     ) 

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
  i = ifelse(value[[1]] > 10, "#E41A1C", ifelse(value[[1]] > -10, "#D3D3D3", "#377EB8"))
  for(h in seq(-100, 100, by = 25)) {
    circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
  }
  circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
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

circos.par("track.height" = 0.15)

circos.genomicDensity(bed_list, col = c("#FF000080", "#0000FF80"))

circos.genomicLabels(annot, labels.column=5, cex=0.5, col="black", line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.05, labels_height=0.04)

dev.off()



###############################

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
                      i = ifelse(value[[1]] > 10, "#E41A1C", ifelse(value[[1]] > -10, "#D3D3D3", "#377EB8"))
                      for(h in seq(-100, 100, by = 25)) {
                        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#E7E7E7")
                      }
                      circos.genomicPoints(region, value, cex = cex, pch = 16, col = i, ...)
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

circos.par("track.height" = 0.15)

circos.genomicDensity(bed_list, col = c("#FF000080", "#0000FF80"))

circos.genomicLabels(annot, labels.column=5, cex=0.5, col= annot$color, line_lwd=0.5, line_col="grey80", 
                     side="inside", connection_height=0.05, labels_height=0.04)
