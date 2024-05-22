if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("cowplot", quietly = TRUE)){
  install.packages("cowplot")
}

if (!require("ggplotify", quietly = TRUE)){
  install.packages("ggplotify")
}

if (!require("ggpubr", quietly = TRUE)){
  install.packages("ggpubr")
}

if (!require("extrafont", quietly = TRUE)){
  install.packages("extrafont")
}

if (!require("methylKit", quietly = TRUE)){
  BiocManager::install("methylKit")
}

if (!require("clusterProfiler", quietly = TRUE)){
  BiocManager::install("clusterProfiler")
}

# circular plot


meth_diff <- readRDS("./data/meth_diff_cut_10_tiles_100.rds")
qvalue_cut <- 0.05
meth_cut <- 10
meth_diff_cut_hypo <- getMethylDiff(meth_diff,
                                    difference = meth_cut,
                                    qvalue = qvalue_cut,
                                    type = "hypo")

meth_diff_cut_hyper <- getMethylDiff(meth_diff,
                                     difference = meth_cut,
                                     qvalue = qvalue_cut,
                                     type = "hyper")
feature_sum_annot <- readRDS("./data/annotated_df__cut_10_tiles100.rds")
file_path_name <- "cut_10_tiles_100"

source("code/circular_plot.R")


circular_plot <- function() {
  circos.clear()
  
  circos.par(gap.after = c(rep(1.5, 9), 12), start.degree = 90)
  
  circos.initializeWithIdeogram(plotType = c("labels", "axis"),
                                species = "xenTro10"
  )
  circos.par("track.height" = 0.5)
  
  
  circos.genomicTrack(DMR_all2[,1:4], 
                      stack = FALSE, 
                      ylim = c(-100, 100),
                      panel.fun = function(region, value, ...) {
                        # cex = DMR_all2[rownames(value), "cex"]
                        # i = DMR_all2[rownames(value), "color"]
                        cex = 0.5
                        i = DMR_all[rownames(value), "color"]
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
  # circos.genomicLabels(annot, labels.column = 5, cex = 1.2, col= annot$pointcolor, line_lwd = 1, line_col="grey10", 
  #                      side="inside", connection_height=0.025, labels_height=0.05, padding = 0.8)
  # 
}

svglite(paste0("./figures/circular_plot_clean.svg"),
        width = 10,
        height = 10,
        scaling = 2)

circular_plot()

dev.off()

# genomic locations plot

meth_diff_cut_annotated <- readRDS("data/meth_diff_cut_annotated_cut_10_tiles100.rds")
plot_genomic_location <- function() {genomation::plotTargetAnnotation(meth_diff_cut_annotated, main = "Genomic Features")}

gen_regions <- genomation::getTargetAnnotationStats(meth_diff_cut_annotated, percentage=FALSE,precedence=TRUE)

reg_perc <- meth_diff_cut_annotated@precedence %>% 
  data.frame %>%
  rownames_to_column(., var = "location")
reg_perc$location <- 
  factor(reg_perc$location, 
         levels = c("promoter", "exon", "intron", "intergenic"))

pie_df <- reg_perc %>% 
  arrange(desc(location)) %>%
  mutate(prop = . / sum(reg_perc$.) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(round = round(., digits = 0))

pie_graph <- ggplot(pie_df, aes(x="", y=prop, fill=location)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="right") +
  # geom_text(aes(y = ypos, 
  #               label = paste(round, "%")), 
  #           color = "black", size=5,
  #           family = "sans"
  #           ) +
  geom_label(aes(y = ypos, label = paste(round, "%")),
             alpha = 1,
             label.padding = unit(0.25, "lines"),
             label.r = unit(0.15, "lines"),
             label.size = 0.25,
             nudge_x = 0.5,
             show.legend = FALSE,
             fill = "white"
  ) +
  scale_fill_brewer(palette="Dark2"
                    ) +
  theme(legend.title=element_blank(),
        legend.position = "bottom")


# convert base graph to grob

c_grob <- as_grob(circular_plot)
gloc_grob <- as_grob(plot_genomic_location)

c_grob <- as.ggplot(circular_plot)

c_grob <- c_grob +
  annotation_custom(grid::textGrob(paste0((nrow(meth_diff_cut_hyper)+nrow(meth_diff_cut_hypo)),
                                          " DMRs \n ", 
                                          nrow(meth_diff_cut_hyper), 
                                          " hypermethylated \n ",
                                          nrow(meth_diff_cut_hypo),
                                          " hypomethylated" ),
                                   gp = grid::gpar(fontsize = 10)), 
                    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) 

# cowplot

figure1 <- plot_grid(c_grob, pie_graph, labels = c("A","B"),
                     scale = c(1.1,0.8)
                     )

ggsave2("figures/Figure_1.tiff",
        plot = figure1,
        scale = 0.8,
        width = 30,
        height = 15,
        units = "cm",
        dpi = 300,
        bg = "white"
)

loadfonts()


ggsave2("figures/Figure_1.eps",
        plot = figure1,
        device = cairo_ps,
        scale = 0.8,
        width = 30,
        height = 15,
        units = "cm",
        dpi = 300,
        bg = "white"
)


##################
# Over-representation plot

go_hypo <- readRDS("data/go_overrep_hypo__cut_10_tiles100.rds")
p_go_hypo <- barplot(go_hypo, 
                     color = "qvalue", 
                     showCategory = 10, 
                     title = paste0("GO Over-representation: \nhypomethylated regions")) +
  xlim(0,7) +
  xlab(NULL) +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      limits=c(0,0.05), 
                      labels=c(0,0.025,0.05),
                      breaks=c(0,0.025,0.05)) 
  # scale_y_discrete(labels = scales::label_wrap(width = 40))

go_prom_hypo <- readRDS("data/go_overrep_prom_hypo_cut_10_tiles_100.rds")
p_go_prom_hypo <- barplot(go_prom_hypo, 
                          color = "qvalue", 
                          showCategory = 10, 
                          title = paste0("GO Over-representation: \nhypomethylated promoters")) +
  xlim(0,7) +
  xlab(NULL) +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      limits=c(0,0.05), 
                      labels=c(0,0.025,0.05),
                      breaks=c(0,0.025,0.05)) 
  # scale_y_discrete(labels = scales::label_wrap(width = 40))



kegg_hyper <- readRDS("data/kegg_overrep_hyper_cut_10_tiles_100.rds")

p_kegg_hyper <- kegg_hyper %>% 
  filter(., Count > 0) %>% 
  arrange(desc(Count)) %>% 
  barplot(., showCategory = 20, 
          title = paste0("KEGG Pathways Over-representation: \nhypermethylated regions")) +
  xlim(0,7) +
  scale_fill_gradient(low = "blue", 
                      high = "red", 
                      limits=c(0,0.05), 
                      labels=c(0,0.025,0.05),
                      breaks=c(0,0.025,0.05)) 
  # scale_y_discrete(labels = scales::label_wrap(width = 30)) 


figure2 <- ggarrange(p_go_hypo, p_go_prom_hypo, p_kegg_hyper, 
          ncol = 1, 
          nrow = 3, 
          heights = c(0.25, 0.65, 0.2),
          common.legend = TRUE, 
          legend = "right",
          labels = "AUTO",
          align = "hv")


ggsave2("figures/Figure_2.tiff",
        plot = figure2,
        scale = 1,
        width = 15,
        height = 25,
        units = "cm",
        dpi = 300,
        bg = "white"
)

ggsave2("figures/Figure_2.eps",
        plot = figure2,
        device = cairo_ps,
        scale = 1,
        width = 15,
        height = 25,
        units = "cm",
        dpi = 300,
        bg = "white"
)


