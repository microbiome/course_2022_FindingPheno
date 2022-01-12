################################################################################
#####**********************************************************************#####
#####********************** INSTALLATION SCRIPT ***************************##### 
#####**********************************************************************#####
################################################################################

if (!requireNamespace("BiocManager")) {
    install.packages("BiocManager")
}

# List of packages that we need from cran and bioc 
packages <- c("ggplot2", "pheatmap", "stringr", "igraph", "ANCOMBC",
              "microbiome", "httpuv", "microbiomeDataSets", "mia", "caret", "ranger",
              "dplyr", "miaViz", "knitr", "kableExtra", "vegan", "ecodist", "biclust",
              "patchwork", "pdp", "MLmetrics", "precrec")

# The following tries to load all required packages, and if they are not available, installs them.
new.pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]

if (length(new.pkgs)) {
    BiocManager::install(new.pkgs, ask=T)
}

# Print check
sapply(packages, require, character.only = TRUE)