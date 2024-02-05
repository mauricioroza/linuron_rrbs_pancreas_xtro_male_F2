## Filtering

# min_n <- 4L #minimum number of samples per group

# Filter step = Discards bases with coverage below 10 reads and above 99.9th percentile
myobj_filtered <- filterByCoverage(myobj,
                                     lo.count=10,
                                     lo.perc=NULL,
                                     hi.count=NULL,
                                     hi.perc=99.9)
# Normalization
myobj_filtered_norm <- normalizeCoverage(myobj_filtered, method = "median")

# Merge data
myobj_united <- methylKit::unite(myobj_filtered_norm, destrand=FALSE, min.per.group = min_n) #min.per.group = minimum number of samples per replicate needed to cover a region/base

meth_united_dir <- paste0("./data/meth_united_",file_path_name,".rds")
saveRDS(myobj_united, ascii=FALSE, file = meth_united_dir)

# Number of regions covered
n_cpg <- nrow(myobj_united)

##Further Filtering

# get percent methylation matrix
perc_meth_matrix <- percMethylation(myobj_united)

# calculate standard deviation of CpGs
sdev_meth <- matrixStats::rowSds(perc_meth_matrix, na.rm = TRUE)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sdev_meth, breaks = 100, xlab = "Methylation Standard Deviation (%)", main = "Histogram of % Methylation Standard Deviation")

# keep only CpG with standard deviations larger than 2%
meth_filtered <- myobj_united[sdev_meth > 2]
meth_filtered <- meth_filtered[!grepl("Sca", meth_filtered$chr), ] #remove CpG located in scaffolds, to smooth further analysis

n_cpg_filtered <- nrow(meth_filtered)

# remove big objects to free memory
# rm(list = c("myobj", "myobj_filtered", "myobj_filtered_norm", "myobj_united"))