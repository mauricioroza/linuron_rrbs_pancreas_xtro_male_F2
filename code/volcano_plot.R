
plot_table <- merge(getData(meth_diff) , feature_sum_annot, all.x = TRUE) %>%
  mutate(methylation_status = case_when(
    meth.diff > meth_cut & qvalue < qvalue_cut ~ "Hypermethylated",
    meth.diff < -meth_cut & qvalue < qvalue_cut ~ "Hypomethylated",
    TRUE ~ "N.S."
  ))

#ggplot
volcano_plot <- ggplot(data=plot_table, aes(x=meth.diff, y=-log10(qvalue), col= methylation_status, fill = methylation_status)) +
  geom_point(shape = 21, size = 2, stroke = 0.1) + 
  theme_minimal() +
  scale_color_manual(values=c("black", "black", "black")) +
  scale_fill_manual(values=c("#E41A1C","#377EB8", "gray")) +
  geom_vline(xintercept=c(-meth_cut, meth_cut), col="gray57", linetype = "longdash") +
  geom_hline(yintercept=-log10(qvalue_cut), col="gray57", linetype = "longdash") +
  labs(title = analysis_name,
       subtitle =paste0("qvalue < ",
                        as.character(qvalue_cut),
                        " & methylation difference >= ",
                        "±",
                        as.character(meth_cut), "%")) +
  xlab("% Methylation Change")

volcano_plot

ggsave(filename = paste0("volcano_plot_", file_path_name, ".tiff"),
       path = "./figures/",
       plot = volcano_plot,
       device = "tiff",
       scale = 2, 
       dpi = 300,
       width = 15,
       height = 10,
       units = "cm",
       bg = "white")
