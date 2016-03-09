setwd("C:/Github_folder/scVEGs")
source("scVEGs.r")
data <- read.delim('data_GSE65525.txt', header = TRUE, stringsAsFactors = TRUE)
rownames(data) <- data[, 1]
data <- data[, -1]
pVal <- 0.01
pFlag <- 1
species <- 'mmu'

cellSum <- apply(data, 2, sum)
scaleNum <- mean(cellSum)/cellSum
scaleFactor <- t(kronecker( matrix(1, 1, dim(data)[1]), scaleNum))
normData <- scaleFactor * data
outputName <- 'result_GSE65525'
sig <- scVEGs(normData, pVal, pFlag, species, outputName)

write.table(sig, file = 'result_GSE65525_sig.txt', append = FALSE, quote = FALSE, 
            sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))

