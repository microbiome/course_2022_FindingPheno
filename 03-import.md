# Data

This section demonstrates how to import data in R.

## Data structure

Such analysis using the miaverse framework, are based upon core data structures 
including SingleCellExperiment (SCE), SummarizedExperiment (SE), TreeSummarizedExperiment (TreeSE) and MultiAssayExperiment (MAE) [(resources)](https://microbiome.github.io/course_2022_miaverse/study-material.html#resources-for-treesummarizedexperiment).

Multi-assay data can be stored in [altExp](https://microbiome.github.io/OMA/containers.html#alternative-experiments) 
slot of TreeSE or [MAE](https://microbiome.github.io/OMA/containers.html#multiassayexperiments) data container.

Different data sets are first imported into SE or TreeSE data container similarly to the case when only one data set is present. After that different data sets are combined into the same data container. Result is one TreeSE object with alternative experiment in altExp slot, or MAE object with multiple experiment in its experiment slot.

## Example data

As an example data, we use data from following publication: Hintikka L et al. (2021) [Xylo-oligosaccharides in prevention of hepatic steatosis and adipose tissue inflammation: associating taxonomic and metabolomic patterns in fecal microbiotas with biclustering](https://www.mdpi.com/1660-4601/18/8/4049).

This example data can be loaded from [microbiomeDataSets](https://bioconductor.org/packages/release/data/experiment/html/microbiomeDataSets.html). The data is already in MAE format. It includes three different experiments: microbial abundance data, metabolite concentrations, and data about different biomarkers.

## Importing data in R


```r
library(stringr)

# Load the data
mae <- microbiomeDataSets::HintikkaXOData()

# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]

# Clean taxonomy data, so that names do not include addtional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))

mae
```

```
## A MultiAssayExperiment object of 3 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 3:
##  [1] microbiota: SummarizedExperiment with 12613 rows and 40 columns
##  [2] metabolites: SummarizedExperiment with 38 rows and 40 columns
##  [3] biomarkers: SummarizedExperiment with 39 rows and 40 columns
## Functionality:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DataFrame
##  sampleMap() - the sample coordination DataFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DataFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
##  exportClass() - save data to flat files
```
