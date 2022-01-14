

# Cross-correlation & Unsupervised learning

## Cross-correlation

With cross-correlation analysis, we can analyze how strongly and how differently variables 
are associated between each other. For instance, we can analyze if higher presence of a 
specific taxon equals to higher levels of a biomolecule.

**TASK**

1. Run installation script to load packages into the session

2. Import HintikkaXO data

3. Subset _microbiota_ data (rank = "Phylum", prevalence = 0.2, detection = 0.001) ([subsetByPrevalentTaxa](https://microbiome.github.io/OMA/differential-abundance.html#prevalence-filtering))

4. Apply clr-transform to _microbiota_ data and log10-tranform to _metabolites_ data ([transformSamples](https://microbiome.github.io/OMA/taxonomic-information.html#data-transformation))

5. Remove uncultured and ambiguous taxa (as it's hard to interpret their results) ( USE THIS: _mae[[1]] <- mae[[1]][-grep("uncultured|Ambiguous_taxa", names(mae[[1]])),]_ )

5. Calculate cross-correlation between _microbiota_ (clr) and _metabolites_ (log10) (Use _show_warnings = FALSE_ as an argument)
([testExperimentCrossCorrelation](https://microbiome.github.io/OMA/multi-assay_analyses.html#multi-assay_analyses))

6. Create a heatmap from cross-correlation matrix [pheatmap](https://microbiome.github.io/OMA/microbiome-community.html#composition-heatmap)

















## Unsupervised learning
Unsupervised learning is a part of machine learning where we try to find information
from unknown data. It is also called data mining. Usually this means finding of
clusters, for instance. Cluster refers to group of samples/features that are similar
between each other. For example, based on clinical data we can try to find patient 
groups that have similar response to used drug.

### Biclustering

Biclustering is a clustering method, which simultaneously clusters rows and columns.
In this example, the aim is to find clusters where subset of taxa share similar
pattern over subset of metabolites. In our case, we try to find clusters where 
taxa and metabolites correlate similarly. 

Check more from OMA which has dedicated chapter on 
[biclustering](https://microbiome.github.io/OMA/biclustering.html). 



```r
# Load package
library(biclust)
```

```
## Loading required package: MASS
```

```
## Loading required package: grid
```

```
## 
## Attaching package: 'grid'
```

```
## The following object is masked from 'package:Biostrings':
## 
##     pattern
```

```
## Loading required package: colorspace
```

```
## Loading required package: lattice
```

```r
# Find biclusters
bc <- biclust(corr$cor, method=BCPlaid(), fit.model = y ~ m)
```

```
## layer: 0 
##  0.07131276
## layer: 1 
## [1]  0 10 15
## [1]  1  8 12
## [1]  2  8 10
## [1] 3 8 8
## [1] 30  8  8
## [1] 31  5  8
## [1] 32  5  6
## [1] 33  4  6
## [1] 34  4  6
## [1] 35  4  6
## [1] 60  4  6
## [1] 4
## [1] 3.75515 0.00000 0.00000 0.00000
## back fitting 2 times
## layer: 2 
## [1]  0 11 17
## [1]  1  8 13
## [1]  2  7 12
## [1] 30  7 12
## [1] 31  0 12
## [1] 32
## [1] 0 0 0 0
##      
## Layer Rows Cols Df   SS   MS Convergence Rows Released Cols Released
##     0   23   38  1 0.00 0.00          NA            NA            NA
##     1    4    6  1 3.97 3.97           1             4             2
```

```r
bc
```

```
## 
## An object of class Biclust 
## 
## call:
## 	biclust(x = corr$cor, method = BCPlaid(), fit.model = y ~ m)
## 
## There was one cluster found with
##   4 Rows and  6 columns
```


```r
# Get biclusters
bicluster_rows <- bc@RowxNumber
bicluster_columns <- bc@NumberxCol

# Convert them into data.frames
bicluster_rows <- data.frame(bicluster_rows)
bicluster_columns <- data.frame(t(bc@NumberxCol))

# Adjust names of clusters
colnames(bicluster_rows) <- paste0("cluster_", 1:ncol(bicluster_rows))
colnames(bicluster_columns) <- paste0("cluster_", 1:ncol(bicluster_columns))

# Print biclusters for rows
head(bicluster_rows)
```

```
##   cluster_1
## 1     FALSE
## 2     FALSE
## 3     FALSE
## 4     FALSE
## 5     FALSE
## 6     FALSE
```

Now, we can add bicluster information into the heatmap that we already made.


```r
# Convert boolean values into numeric
bicluster_columns[ , 1] <- as.numeric(bicluster_columns[, 1] )
bicluster_rows[ , 1] <- as.numeric(bicluster_rows[ , 1])

# Adjust their rownames
rownames(bicluster_columns) <- colnames(corr$cor)
rownames(bicluster_rows) <- rownames(corr$cor)

# Get correlation values that are over thresholds
p_threshold <- 0.01
corr_values <- ifelse(corr$p_adj<p_threshold, round(corr$cor,1), "")

# Create a heatmap
pheatmap(corr$cor,
         annotation_col = bicluster_columns, 
         annotation_row = bicluster_rows,
         
         display_numbers = corr_values,
   
         main = paste0("Correlations between bacteria and metabolites 
                       (correlation over threshold (p < ", p_threshold,") marked)"),
         
         breaks = breaks,
         color = colors, 
   
         fontsize_number = 6, 
         number_color = "yellow")
```

![](07-unsupML_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 
