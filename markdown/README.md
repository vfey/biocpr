### About BioCPR

Correlation analysis is a powerful approach to determine the strength of association between two variables. The strength and direction of a correlation is represented in the form of a correlation coefficient ranging from -1 (negative) to +1 (positive). In cases of biological data such as expression data, establishing a clear pattern can be difficult from just looking at the values. The BioCPR tool has been created to enable researchers to see the results of correlation analysis graphically in the form of an interactive heatmap, which can aid in better understanding of relationships and help identify new patterns.

### Disclaimer

BioCPR tool is intended for visualization purpose only and should not be used for treating or diagnosing human subjects.

### Availability  

The source code and installation instructions can be obtained from [gitlab](https://github.com/vfey/biocpr/).

### Required Software

BioCPR requires 4 R-packages in-order to run. 

You can either install the stable version from CRAN and it will install all the necessary depencies for running the tool. Please find the commands below;

install.packages("heatmapFlex")
install.packages("convertid")
install.packages("readmoRe")
install.packages("coreheat")

If you prefer to install the latest development version, there are 4 custom in-built libraries that needs to be installed in addition to libraries from CRAN and BioConductor. Please find the commands below;

.libPaths()
libLocation <- .libPaths()[1]

## CRAN packages ##
cranPackages <- c("shinyBS", "shinythemes", "knitr", "rmarkdown", "shinyjs", "plyr", "RColorBrewer", "R.utils", "gdata", "data.table", "foreach", "ggplot2", "scales", "curl", "openssl", "httr", "Rcurl", "XML", "WGCNA", "DT", "devtools")
sapply(cranPackages, install.packages, lib = libLocation)

## Bioconductor packages ##

bioConPackages <- c("AnnotationDbi", "biomaRt", "org.Hs.eg.db", "org.Mm.eg.db", "Heatplus", "genefilter", "impute", "preprocessCore", "GO.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(bioConPackages)

## In-house packages ##

install_github("vfey/heatmapFlex")
install_github("vfey/convertid")
install_github("vfey/readmoRe")
install_github("vfey/coreheat")

The required packages can also be obtained from this [link](https://github.com/vfey/biocpr/tree/main/packages/).

### Frequently asked questions  

Please refer to the [FAQ](https://github.com/vfey/biocpr/blob/main/markdown/FAQ.md) here.
