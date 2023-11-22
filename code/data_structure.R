
##Data Structure/Outlier Detection

# Correlation matrix
# getCorrelation(meth_filtered, plot = FALSE)

meth_corr <- cor(percMethylation(meth_filtered), use = "complete.obs", method = "pearson")

meth_corr_longer <- meth_corr %>% data.frame %>% 
  rownames_to_column(var = "rowname") %>% 
  pivot_longer(-rowname, names_to = "colname", values_to = "value") %>%
  mutate(colname = str_remove_all(colname, "X"))
  
corrplot <- ggplot(data = meth_corr_longer, aes(x = rowname, y = colname, fill = value)) + 
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's r")

# Hierarchical clustering
cluster_meth <- function() {clusterSamples(meth_filtered, dist = "correlation", method="ward.D", plot=TRUE)}
cluster_meth()

# Principal Component Analysis
PCA_meth <- function() {PCASamples(meth_filtered, comp = c(1,2))}
PCA_meth()