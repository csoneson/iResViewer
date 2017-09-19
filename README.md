## iResViewer

 <span style="color:blue">**Note: This is work in development, mainly for personal use, and it is provided "as-is", without guarantees of robustness, correctness or optimality. Please use responsibly.**</span>

`iResViewer` is an R package for generating an interactive result visualization for a differential gene expression analysis. 

### Installation

```
## List dependencies
pkg <- c("shiny", "ggplot2", "ggrepel", "dplyr", "tidyr", "Gviz", 
         "DT", "shinydashboard", "GenomicRanges", "IRanges", 
         "rtracklayer", "S4Vectors", "GenomeInfoDb", 
         "BiocGenerics", "devtools")

## Check if dependencies are already installed
pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]

## If some dependency is missing, install it
if (length(pkg) > 0) {
	source("https://bioconductor.org/biocLite.R")
	biocLite(pkg, dependencies = TRUE, ask = FALSE)
}

## Install iResViewer
devtools::install_github("csoneson/iResViewer")
```