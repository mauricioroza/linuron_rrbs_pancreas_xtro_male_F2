# check number of samples
lapply(myobj, function(x) {nrow(x)})

# What type of data is stored here?
head(myobj[[1]])

# Get a histogram of the methylation percentage per sample
par(mfrow= c(3,3))
lapply(myobj, function(x) {getMethylationStats(x, plot=TRUE, both.strands=FALSE)})

# Get a histogram of the read coverage per sample
par(mfrow= c(3,3))
lapply(myobj, function(x) {getCoverageStats(x, plot=TRUE, both.strands=FALSE)})
