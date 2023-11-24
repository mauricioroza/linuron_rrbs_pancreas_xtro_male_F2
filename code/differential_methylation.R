if (!require("DT", quietly = TRUE)){
  install.packages("DT")
}

# Test for differential methylation... This might take a few minutes.
meth_diff <- calculateDiffMeth(meth_filtered,
                               overdispersion = "MN",
                               adjust="SLIM",
                               test = "Chisq"
                               )


#Save data in RData file
meth_diff_dir <- paste0("./data/meth_diff_",file_path_name,".rds")
saveRDS(meth_diff, ascii=FALSE, file = meth_diff_dir)

# meth_diff <- meth_diff[!grepl("Sca", meth_diff$chr), ] #remove DMR located in scaffolds of the genome, to smooth further analysis
nrow_meth_diff <- nrow(meth_diff)

#Histograms of pvalue distribution
# hist(meth_diff$pvalue, main = "p-value Distribution")
# hist(meth_diff$qvalue, main = "q-value Distribution")
# hist(meth_diff$meth.diff, main = "% Methylation Distribution")

# cutoff values of % methylation and qvalue
# meth_cut <- 10
# qvalue_cut <- 0.05
# change directly on source("code/methylation_qvalue_cut_values.R")

# get all differentially methylated bases and order by qvalue
meth_diff_cut <- getMethylDiff(meth_diff,
                               difference = meth_cut,
                               qvalue = qvalue_cut)

meth_diff_cut_hypo <- getMethylDiff(meth_diff,
                               difference = meth_cut,
                               qvalue = qvalue_cut,
                               type = "hypo")

meth_diff_cut_hyper <- getMethylDiff(meth_diff,
                                    difference = meth_cut,
                                    qvalue = qvalue_cut,
                                    type = "hyper")

meth_diff_cut <- meth_diff_cut[order(meth_diff_cut$qvalue),]

nrow_meth_diff_cut <- nrow(meth_diff_cut)

nrow_hypo <- nrow(meth_diff_cut_hypo)

nrow_hyper <- nrow(meth_diff_cut_hyper)
