################################################################################
#####**********************************************************************#####
#####*********************** DAY 2 PRACTICAL ******************************##### 
#####**********************************************************************#####
################################################################################

################################################################################
# Loading data and pre-processing
################################################################################

library(mia)
library(miaViz)
library(dplyr)
library(stringr)
mae <- microbiomeDataSets::HintikkaXOData()
# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]
# Clean taxonomy data, so that names do not include addtional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))
# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                                colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")
# Calculates relative abundances, and stores the table to assays
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")
# agglomerate the microbiota experiment to Phylum
se_phylum <- agglomerateByRank(mae[[1]], rank = "Phylum")


################################################################################
# Chapter 6 Beta diversity
################################################################################

### 6.1.1 PCoA for ASV-level data with Bray-Curtis
# Pick the relative abundance table
rel_abund_assay <- assays(mae[[1]])$relabundance

# Calculates Bray-Curtis distances between samples. Because taxa is in
# columns, it is used to compare different samples. We transpose the
# assay to get taxa to columns
bray_curtis_dist <- vegan::vegdist(t(rel_abund_assay), method = "bray")

# PCoA
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

# All components could be found here: 
# bray_curtis_pcoa$vectors
# But we only need the first two to demonstrate what we can do:
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])

# Create a plot
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
    geom_point() +
    labs(x = "PC1",
         y = "PC2", 
         title = "Bray-Curtis PCoA") +
    theme_bw(12) # makes titles smaller

bray_curtis_plot

## 6.2 Highlighting external variables

# Add diet information to data.frame
bray_curtis_pcoa_df$Diet <- colData(mae)$Diet

# Creates a plot
plot <- ggplot(data = bray_curtis_pcoa_df, aes_string(x = "pcoa1", y = "pcoa2", color = "Diet")) +
    geom_point() +
    labs(x = "PC1",
         y = "PC2", 
         title = "Bray-Curtis PCoA") +
    theme_bw(12) 

plot

## 6.3 Estimating associations with an external variable

# First we get the relative abundance table
rel_abund_assay <- assays(mae[[1]])$relabundance

# again transpose it to get taxa to columns
rel_abund_assay <- t(rel_abund_assay)

# then we can perform the method
permanova_diet <- vegan::adonis(rel_abund_assay ~ Diet,
                                data = colData(mae),
                                permutations = 99)

# we can obtain a the p value for our predictor:
print(paste0("The test result p-value: ", 
             as.data.frame(permanova_diet$aov.tab)["Diet", "Pr(>F)"]))

# Gets the coefficients
coef <- coefficients(permanova_diet)["Diet1",]

# Gets the highest coefficients
top.coef <- sort(head(coef[rev(order(abs(coef)))],20))

# Plots the coefficients
top_taxa_coeffient_plot <- ggplot(data.frame(x = top.coef,
                                             y = factor(names(top.coef),
                                                        unique(names(top.coef)))),
                                  aes(x = x, y = y)) +
    geom_bar(stat="identity") +
    labs(x="", y="", title="Top Taxa") +
    theme_bw()

top_taxa_coeffient_plot

# Gets corresponding Genus level names and stores them to top.coef
names <- rowData(mae[[1]])[names(top.coef), ][,"Genus"]

# Adds new labels to the plot
top_taxa_coeffient_plot <- top_taxa_coeffient_plot +
    scale_y_discrete(labels = names) # Adds new labels
top_taxa_coeffient_plot

bray_curtis_pcoa_df$Diet <- colData(mae)$Diet
p_values <- list()
for(pc in c("pcoa1", "pcoa2")){
    # Creates a formula from objects
    formula <- as.formula(paste0(pc, " ~ ", "Diet"))
    # Does the permanova analysis
    p_values[[pc]] <- vegan::adonis(formula, data = bray_curtis_pcoa_df,
                                    permutations = 9999, method = "euclidean"
    )$aov.tab["Diet", "Pr(>F)"]
}

# Creates a plot
plot <- ggplot(data = bray_curtis_pcoa_df, aes_string(x = "pcoa1", y = "pcoa2", color = "Diet")) +
    geom_point(size = 3) +
    labs(title = paste0("PCoA beta diversity ordination for microbiome samples"),
         x = paste0("PC1 (p = ", p_values[["pcoa1"]], ")"), y = paste0("PC2 (p = ", p_values[["pcoa2"]], ")")) +
    theme_bw() 

plot

################################################################################
# Chapter 7 Unsupervised learning
################################################################################

## 7.1 Biclustering

library(ggplot2)
# Threshold: metabolites whose (cv > +threshold or cv < -threshold), will be included
cv_threshold <- 0.5
metabolite_trans <- "nmr"

# Get the data
metabolite_tse <- mae[[2]]

# Calculate coefficient of variation of individual metabolites
df <- data.frame(cv = apply(assay(metabolite_tse, metabolite_trans), 1, 
                            function(x){sd(x)/mean(x)}))

# Plot them as a histogram, and show a line that is used as a threshold
plot <- ggplot(df, aes(x = cv)) +
    geom_histogram(bins = 50, color="darkred", fill="lightblue") +
    labs(x = "CV", y = "metabolite frequency", 
         title = "Distribution of coefficient of 
       variation of log10 concentration of metabolites") +
    geom_vline(xintercept = cv_threshold, color = "red") +
    geom_text(aes(cv_threshold, 6, label = 
                      paste0("CV threshold (", cv_threshold, ")"), vjust = 2, angle=90)) +
    geom_vline(xintercept = -cv_threshold, color = "red") +
    geom_text(aes(-cv_threshold, 6, label = 
                      paste0("CV threshold (", -cv_threshold, ")"), vjust = -1, angle=90))

plot

# Get those metabolites that are over threshold
metabolites_over_th <- rownames(df[df$cv > cv_threshold | 
                                       df$cv < -cv_threshold, , drop = FALSE])
# Ignore those metabolites that do not have name / are NA
metabolites_over_th <- metabolites_over_th[!str_detect(metabolites_over_th, "NA")]

rank <- "Genus"
prevalence <- 0.2
detection <- 0.001
taxa_trans <-  "clr"

# Get bacterial data
taxa_tse <- mae[[1]]
# Agglomerate at Genus level
taxa_tse <- agglomerateByRank(taxa_tse, rank = rank)
# Do CLR transformation
taxa_tse <- transformSamples(taxa_tse, method = "clr", pseudocount = 1)

# Subset metabolite data
metabolite_tse <- metabolite_tse[metabolites_over_th, ]

# Subset bacterial data by its prevalence. Bacteria whose prevalences are over 
# threshold are included
taxa_tse <- subsetByPrevalentTaxa(taxa_tse, 
                                  prevalence = prevalence, 
                                  detection = detection)

# Remove uncultured and ambiguous(as it's hard to interpret their results)
taxa_tse <- taxa_tse[-grep("uncultured|Ambiguous_taxa", names(taxa_tse)),]

library(pheatmap)

# Cross correlate data sets
correlations <- testExperimentCrossCorrelation(taxa_tse, metabolite_tse, 
                                               abund_values1 = "clr", abund_values2 = "nmr",
                                               method = "spearman", mode = "matrix")

# For plotting purpose, convert p-values, under 0.05 are marked with "X"
p_threshold <- 0.01
p_values <- ifelse(correlations$p_adj<p_threshold, "X", "")

# Scale colors
breaks <- seq(-ceiling(max(abs(correlations$cor))), ceiling(max(abs(correlations$cor))), 
              length.out = ifelse( max(abs(correlations$cor))>5, 
                                   2*ceiling(max(abs(correlations$cor))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", 
                             "red", "darkred"))(length(breaks)-1)

# Create a heatmap
pheatmap(correlations$cor, display_numbers = p_values,
         main = paste0("Correlations between bacteria and metabolites 
              (statistically significant associations (p < 0.05) marked with X)"),
         fontsize = 10,
         breaks = breaks,
         color = colors, 
         fontsize_number = 10)

# Load package
library(biclust)

# Find biclusters
bc <- biclust(correlations$cor, method=BCPlaid(), fit.model = y ~ m,
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10,
              iter.startup = 10, iter.layer = 100, verbose = FALSE)

bc

# Functions for obtaining biclust information

# Get clusters for rows and columns
.get_biclusters_from_biclust <- function(bc, assay){
    # Get cluster information for columns and rows
    bc_columns <- t(bc@NumberxCol)
    bc_columns <- data.frame(bc_columns)
    bc_rows <- bc@RowxNumber
    bc_rows <- data.frame(bc_rows)
    
    # Get data into right format
    bc_columns <- .manipulate_bc_data(bc_columns, assay, "col")
    bc_rows <- .manipulate_bc_data(bc_rows, assay, "row")
    
    return(list(bc_columns = bc_columns, bc_rows = bc_rows))
}

# Input clusters, and how many observations there should be, i.e., the number of samples or features
.manipulate_bc_data <- function(bc_clusters, assay, row_col){
    # Get right dimension
    dim <- ifelse(row_col == "col", ncol(assay), nrow(assay))
    # Get column/row names
    if( row_col == "col" ){
        names <- colnames(assay)
    } else{
        names <- rownames(assay)
    }
    
    # If no clusters were found, create one. Otherwise create additional cluster which
    # contain those samples that are not included in clusters that were found.
    if( nrow(bc_clusters) != dim ){
        bc_clusters <- data.frame(cluster = rep(TRUE, dim))
    } else {
        # Create additional cluster that includes those samples/features that
        # are not included in other clusters.
        vec <- ifelse(rowSums(bc_clusters) > 0, FALSE, TRUE)
        # If additional cluster contains samples, then add it
        if ( any(vec) ){
            bc_clusters <- cbind(bc_clusters, vec)
        }
    }
    # Adjust row and column names
    rownames(bc_clusters) <- names
    colnames(bc_clusters) <- paste0("cluster_", 1:ncol(bc_clusters))
    return(bc_clusters)
}

# Get biclusters
bcs <- .get_biclusters_from_biclust(bc, correlations$cor)

bicluster_rows <- bcs$bc_rows
bicluster_columns <- bcs$bc_columns

# Print biclusters for rows
head(bicluster_rows)

# Convert boolean values into factors
bicluster_columns <- data.frame(apply(bicluster_columns, 2, as.factor))
bicluster_rows <- data.frame(apply(bicluster_rows, 2, as.factor))

# Adjust colors for all clusters
if( ncol(bicluster_rows) > ncol(bicluster_columns) ){
    cluster_names <- colnames(bicluster_rows)
} else {
    cluster_names <- colnames(bicluster_columns)
}
annotation_colors <- list()
for(name in cluster_names){
    annotation_colors[[name]] <- c("TRUE" = "red", "FALSE" = "white")
}

# Get correlation values that are over thresholds
p_threshold <- 0.01
corr_threshold <- 0.6
corr_values <- ifelse(correlations$p_adj<p_threshold & 
                          abs(correlations$cor)>corr_threshold , round(correlations$cor,1), "")

# Create a heatmap
pheatmap(correlations$cor,
         annotation_col = bicluster_columns, 
         annotation_row = bicluster_rows,
         annotation_colors = annotation_colors,
         display_numbers = corr_values,
         main = paste0("Correlations between bacteria and metabolites 
              (correlation over threshold (corr > ", corr_threshold,
                       ", p < ", p_threshold,") marked)"),
         fontsize = 10,
         breaks = breaks,
         color = colors, 
         fontsize_number = 6, 
         number_color = "yellow",
         annotation_legend = FALSE)

################################################################################
# Chapter 8 Supervised learning
################################################################################

# Clearing-up variable space
rm(list = ls())

# Loading data and pre-processing

mae <- microbiomeDataSets::HintikkaXOData()

# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]

# Clean taxonomy data, so that names do not include addtional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2, 
                                     str_remove, pattern = "._[0-9]__"))

# For simplicity, classify all high-fat diets as high-fat, and all the low-fat 
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" | 
                                colData(mae)$Diet == "High-fat + XOS", 
                            "High-fat", "Low-fat")

# Calculates relative abundances, and stores the table to assays
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")

# Earlier processing
# Threshold: metabolites whose (cv > +threshold or cv < -threshold), will be included
cv_threshold <- 0.5
metabolite_trans <- "nmr"

# Get the data
metabolite_tse <- mae[[2]]

# Calculate coeffieicnt of variation of individual metabolites
df <- data.frame(cv = apply(assay(metabolite_tse, metabolite_trans), 1, 
                            function(x){sd(x)/mean(x)}))

# Get those metabolites that are over threshold
metabolites_over_th <- rownames(df[df$cv > cv_threshold | 
                                       df$cv < -cv_threshold, , drop = FALSE])
# Ignore those metabolites that do not have name / are NA
metabolites_over_th <- metabolites_over_th[!str_detect(metabolites_over_th, "NA")]

rank <- "Genus"
prevalence <- 0.2
detection <- 0.001
taxa_trans <-  "rclr"

# Get bacterial data
taxa_tse <- mae[[1]]
# Agglomerate at Genus level
taxa_tse <- agglomerateByRank(taxa_tse, rank = rank)
# Do CLR transformation
taxa_tse <- transformSamples(taxa_tse, method = "rclr", pseudocount = 1)

# Subset metabolite data
metabolite_tse <- metabolite_tse[metabolites_over_th, ]

# Subset bacterial data by its prevalence. Bacteria whose prevalences are over 
# threshold are included
taxa_tse <- subsetByPrevalentTaxa(taxa_tse, 
                                  prevalence = prevalence, 
                                  detection = detection)

# Remove uncultured and ambiguous(as it's hard to interpret their results)
taxa_tse <- taxa_tse[-grep("uncultured|Ambiguous_taxa", names(taxa_tse)),]

# Define data sets to cross-correlate
x <- t(assay(taxa_tse, taxa_trans))
y <- t(assay(metabolite_tse, "nmr"))
# If there are duplicated taxa names, makes them unique
colnames(x) <- str_remove(colnames(x), paste0(rank, ":"))
colnames(x) <- make.unique(colnames(x))

## 8.1 Data curation

butyrate_df <- data.frame(cbind(y, x))
butyrate_df <- butyrate_df[,which(colnames(butyrate_df) %in% c("Butyrate", colnames(x)))]

library(caret)
set.seed(42)
trainIndex <- createDataPartition(butyrate_df$Butyrate, p = .8, list = FALSE, times = 1)
butyrate_df_train <- butyrate_df[trainIndex,]
butyrate_df_test <- butyrate_df[-trainIndex,]

## 8.2 Regression with random forests
set.seed(42)
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 5)
rfFit1 <- train(Butyrate ~ ., data = butyrate_df_train, 
                method = "ranger", 
                trControl = fitControl,
                importance = "permutation")
print(rfFit1)
print(rfFit1$finalModel)

test_predictions <- predict(rfFit1, newdata = butyrate_df_test)
print(postResample(test_predictions, butyrate_df_test$Butyrate))

# Plot predicted vs observed
pred_obs <- data.frame(predicted = test_predictions, observed = butyrate_df_test$Butyrate)
ggplot(data = pred_obs, aes(x=predicted, y=observed)) + geom_point(size = 5, color = "orange") + 
    xlab("Predicted butyrate concentration") + ylab("Observed butyrate concentration") +
    lims(x = c(0,5), y = c(0,5)) +
    geom_abline(linetype = 5, color = "blue", size = 1)+ # Plot a perfect fit line
    theme(panel.border = element_rect(colour = "black", fill = NA),
          panel.background = element_blank())

plot(varImp(rfFit1))

library(patchwork)
library(pdp)
top_features <- rownames(varImp(rfFit1)$importance)[order(varImp(rfFit1)$importance[,"Overall"], decreasing = TRUE)[1:6]]
pd_plots <- list(NULL)
for (feature in 1:length(top_features)) {
    pd_plots[[feature]] <- partial(rfFit1, pred.var = top_features[feature], rug = TRUE) %>% autoplot() + 
        geom_hline(yintercept = mean(butyrate_df_train$Butyrate), linetype = 2, color = "gray") + # Show the mean of the training data as a dashed line
        scale_y_continuous(limits=c(1.5,2.3)) + # Harmonize the scale of yhat on all plots
        theme(panel.border = element_rect(colour = "black", fill = NA),
              panel.background = element_blank())
    print(paste0("Partial dependence of ", top_features[feature]))
}

wrap_plots(pd_plots)

## 8.3 Classification with random forests

butyrate_cutoff <- median(butyrate_df_test$Butyrate)
butyrate_df_test_2 <- butyrate_df_test
butyrate_df_train_2 <- butyrate_df_train
butyrate_df_test_2$Butyrate <- as.factor(ifelse(butyrate_df_test_2$Butyrate >= butyrate_cutoff, "High", "Low"))
butyrate_df_train_2$Butyrate <- as.factor(ifelse(butyrate_df_train_2$Butyrate >= butyrate_cutoff, "High", "Low"))

set.seed(42)
fitControl <- trainControl(method = "repeatedcv", number = 5, repeats = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
rfFit2 <- train(Butyrate ~ ., data = butyrate_df_train_2, 
                method = "ranger", 
                trControl = fitControl,
                importance = "permutation")

test_predictions_2 <- data.frame(obs = butyrate_df_test_2$Butyrate, 
                                 pred = predict(rfFit2, newdata = butyrate_df_test_2), 
                                 predict(rfFit2, newdata = butyrate_df_test_2, type = "prob"))

#Print out the metrics
print(rfFit2)

print(rfFit2$finalModel)

print(twoClassSummary(test_predictions_2, lev = c("High", "Low")))

library(MLmetrics)
library(precrec)
aucs <- evalmod(scores = test_predictions_2$Low, labels = test_predictions_2$obs)
print(aucs)

autoplot(aucs)

plot(varImp(rfFit2))
