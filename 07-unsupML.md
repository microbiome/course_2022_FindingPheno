

# Multi-Omics & Unsupervised learning

## Cross-correlation

With cross-correlation analysis, we can analyze how strongly and how differently variables 
are associated between each other. For instance, we can analyze if higher presence of a 
specific taxon equals to higher levels of a biomolecule.

**TASK**

1. Run installation script to load packages into the session

2. Import HintikkaXO data

3. Subset _microbiota_ data (rank = "Genus", prevalence = 0.2, detection = 0.001) ([subsetByPrevalentTaxa](https://microbiome.github.io/OMA/differential-abundance.html#prevalence-filtering))

4. Apply clr-transform to _microbiota_ data and log10-tranform to _metabolites_ data ([transformSamples](https://microbiome.github.io/OMA/taxonomic-information.html#data-transformation))

5. Remove uncultured and ambiguous taxa (as it's hard to interpret their results) ( USE THIS: _mae[[1]] <- mae[[1]][-grep("uncultured|Ambiguous_taxa", names(mae[[1]])),]_ )

5. Calculate cross-correlation between _microbiota_ (clr) and _metabolites_ (log10) (Use _show_warnings = FALSE_, _test_significance = TRUE_, and _mode = "matrix"_ as an arguments)
([getExperimentCrossCorrelation](https://microbiome.github.io/OMA/multi-assay_analyses.html#multi-assay_analyses))

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
set.seed(19574)

# Find biclusters
bc <- biclust(corr$cor, method=BCPlaid(), fit.model = y ~ m, 
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10, 
              iter.startup = 10, iter.layer = 100, verbose = FALSE)

bc
```

```
## 
## An object of class Biclust 
## 
## call:
## 	biclust(x = corr$cor, method = BCPlaid(), fit.model = y ~ m, 
## 	    background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10, 
## 	    iter.startup = 10, iter.layer = 100, verbose = FALSE)
## 
## There was one cluster found with
##   6 Rows and  7 columns
```


```r
# Get biclusters
bicluster_rows <- bc@RowxNumber
bicluster_columns <- bc@NumberxCol

# Convert into data.frames
bicluster_rows <- as.data.frame(bicluster_rows)
bicluster_columns <- as.data.frame(t(bicluster_columns))

# Adjust names of clusters
colnames(bicluster_rows) <- paste0("cluster_", 1:ncol(bicluster_rows))
colnames(bicluster_columns) <- paste0("cluster_", 1:ncol(bicluster_columns))

# Print biclusters for rows
head(bicluster_rows)
```

```
##   cluster_1
## 1     FALSE
## 2      TRUE
## 3     FALSE
## 4     FALSE
## 5     FALSE
## 6     FALSE
```

Now, we can add bicluster information into the heatmap that we already made.


```r
# Convert boolean values into numeric
bicluster_columns[ , 1] <- as.numeric(bicluster_columns[ , 1])
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
   
         fontsize_number = 8, 
         number_color = "yellow")
```

![](07-unsupML_files/figure-latex/unnamed-chunk-11-1.pdf)<!-- --> 

### Multi-Omics Factor Analysis

We use the R [MOFA2](https://biofam.github.io/MOFA2/installation.html) package
for the analysis, and install corresponding dependencies:


```r
reticulate::install_miniconda(force = TRUE)
```

```
## [1] "/github/home/.local/share/r-miniconda"
```

```r
reticulate::use_miniconda()
reticulate::py_install(packages = c("mofapy2"), pip = TRUE)
```

The `mae` object could be used straight to create the MOFA model:


```r
library(MOFA2)

# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                              colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")

mae[[1]] <- transformCounts(mae[[1]], method = "clr", pseudocount = 1)
assay(mae[[1]], "counts") <- NULL
assay(mae[[2]], "nmr") <- NULL

model <- create_mofa_from_MultiAssayExperiment(mae,
                                               groups = "Diet", 
                                               extract_metadata = TRUE)
model
```

```
## Untrained MOFA model with the following characteristics: 
##  Number of views: 3 
##  Views names: microbiota metabolites biomarkers 
##  Number of features (per view): 12706 38 39 
##  Number of groups: 2 
##  Groups names: High-fat Low-fat 
##  Number of samples (per group): 20 20 
## 
```

Model options could be defined as follows:


```r
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 15
head(model_opts)
```

```
## $likelihoods
##  microbiota metabolites  biomarkers 
##  "gaussian"  "gaussian"  "gaussian" 
## 
## $num_factors
## [1] 15
## 
## $spikeslab_factors
## [1] FALSE
## 
## $spikeslab_weights
## [1] TRUE
## 
## $ard_factors
## [1] TRUE
## 
## $ard_weights
## [1] TRUE
```

Model's training options are defined with the following:


```r
train_opts <- get_default_training_options(model)
head(train_opts)
```

```
## $maxiter
## [1] 1000
## 
## $convergence_mode
## [1] "fast"
## 
## $drop_factor_threshold
## [1] -1
## 
## $verbose
## [1] FALSE
## 
## $startELBO
## [1] 1
## 
## $freqELBO
## [1] 5
```

Preparing and training the model:


```r
model.prepared <- prepare_mofa(
  object = model,
  model_options = model_opts
)
model.trained <- run_mofa(model.prepared)
```

Visualizing the variance explained:


```r
library(patchwork)
wrap_plots(
    plot_variance_explained(model.trained, x="view", y="factor", plot_total = T),
    nrow = 2
) + plot_annotation(title = "Varience Explained per factor and assay",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
```

![](07-unsupML_files/figure-latex/unnamed-chunk-17-1.pdf)<!-- --> 

The top weights for metabolites using the three first factors:


```r
p1 <- plot_top_weights(model.trained,
  view = "metabolites",
  factors = 1:3,
  nfeatures = 10
) + labs(title = "Top weights of the metabolomic assay")

p2 <- plot_top_weights(model.trained,
  view = "microbiota",
  factors = 1:3,
  nfeatures = 10
) + labs(title = "Top weights of the microbiota assay")+
    theme(text = element_text(size = 8))

p1/p2
```

![](07-unsupML_files/figure-latex/unnamed-chunk-18-1.pdf)<!-- --> 

More tutorials and examples of using the package are found at: [link](https://biofam.github.io/MOFA2/tutorials.html)
