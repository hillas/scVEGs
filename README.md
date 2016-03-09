# Introduction 

__scVEGs__ is a novel algorithm for single-cell RNA-seq (scRNA-seq) data to determine significant variably expressed genes (VEGs) using a gene expression variation model (GEVM). It utilizes the relation between coefficient of variation (CV) and average expression level to address the over-dispersion of single-cell data, and its corresponding statistical significance to quantify the variably expressed genes. For more details about the algorithm , read the scVEGs paper that will be published in BMC genomics.

# Installation 

scVEGs is a function in R that can be used by 

```{r,eval=FALSE}
source("scVEGs.r")
```

In the future, scVEGs will be integrated with other tools for scRNA-seq to create an R package.

# Use scVEGs to determine significant variably expressed genes

In the R script __scVEGs_script.r__, it is an thorough example to detect VEGs using the function scVEGs.

First, import a data matrix that contains the UMI read count result of a single-cell RNA-seq data set. We have provided a data set __data_GSE65525__, which is obtained from Gene Expression Omnibus [(GSE65525)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525 ). We also assign some cutoff criterions for detecting VEGs. For the species option, scVEGs calculates the transcripts per million (TPM) of the data and removes genes that have less than 1% of cells that TPM > 1. For now, we have provided the gene length information for __hg19__(hg19_genes_length.tsv) and __mm9_-(mm9_genes_length.tsv). The filtering step will be skipped if species is assigned other than hsa or mmu.

```{r}
data <- read.delim('data_GSE65525.txt', header = TRUE, stringsAsFactors = TRUE)
rownames(data) <- data[, 1]
data <- data[, -1]
pVal <- 0.01		# p-value cutoff
pFlag <- 1			# pFlag: 1, use Benjamini adjusted p-value.  0: use raw p-value.
species <- 'mmu'	# Species of input data, hsa or mmu.
```

The following step is to do a scaling normalization across the data set. You can revise or skip this step as your preference.

```{r}
cellSum <- apply(data, 2, sum)
scaleNum <- mean(cellSum)/cellSum
scaleFactor <- t(kronecker( matrix(1, 1, dim(data)[1]), scaleNum))
normData <- scaleFactor * data
```

At last, we can execute the function scVEGs and export the detected significant VEGs.

```{r}
outputName <- 'result_GSE65525'
sig <- scVEGs(normData, pVal, pFlag, species, outputName)
write.table(sig, file = 'result_GSE65525_sig.txt', append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, qmethod = c("escape", "double"))
```

# Bug reports

Report bugs as issues on our [GitHub repository](https://github.com/hillas/scVEGs/issues) or you can report directly to my email: chenhh@uthscsa.edu.
