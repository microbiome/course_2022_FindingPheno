################################################################################
#####**********************************************************************#####
#####*********************** DAY 2 PRACTICAL ******************************##### 
#####**********************************************************************#####
################################################################################

################################################################################
# Loading data and pre-processing
################################################################################


################################################################################
# Chapter 7 Multi-Omics 6 Unsupervised learning
################################################################################


## ----echo=FALSE, message=FALSE, warning=FALSE-----------------------------------------------------------------------------------
library(mia)
library(miaViz)
library(stringr)
library(pheatmap)


## -------------------------------------------------------------------------------------------------------------------------------
mae <- microbiomeDataSets::HintikkaXOData()
mae


## -------------------------------------------------------------------------------------------------------------------------------
mae[[1]] <- as(mae[[1]], "TreeSummarizedExperiment")
altExp(mae[[1]], "Genus") <- subsetByPrevalentTaxa(mae[[1]], rank = "Genus", prevalence = 0.2, detection = 0.001)
altExp(mae[[1]], "Genus")


## -------------------------------------------------------------------------------------------------------------------------------
altExp(mae[[1]], "Genus") <- transformSamples(altExp(mae[[1]], "Genus"), method = "clr", pseudocount = 1)
mae[[2]] <- transformSamples(mae[[2]], abund_values = "nmr", method = "log10")


## -------------------------------------------------------------------------------------------------------------------------------
altExp(mae[[1]], "Genus") <-altExp(mae[[1]], "Genus")[-grep("uncultured|Ambiguous_taxa", 
                                                            names(altExp(mae[[1]], "Genus"))),]


## -------------------------------------------------------------------------------------------------------------------------------
corr <- getExperimentCrossCorrelation(altExp(mae[[1]], "Genus"), 
                                      mae[[2]], 
                                      "clr", 
                                      "log10", 
                                      show_warnings = FALSE,
                                      test_significance = TRUE,
                                      mode = "matrix")
head(corr$cor,2)


## ----fig.height=10, fig.width=16, message=FALSE, warning=FALSE------------------------------------------------------------------
pheatmap(corr$cor)


## ----fig.height=10, fig.width=14, message=FALSE, warning=FALSE------------------------------------------------------------------
mat <- corr$cor
# Determines the scaling of colors
# Scale colors
breaks <- seq(-ceiling(max(abs(mat))), ceiling(max(abs(mat))), 
              length.out = ifelse( max(abs(mat))>5, 2*ceiling(max(abs(mat))), 10 ) )
colors <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(length(breaks)-1)

# For plotting purpose, convert p-values, under 0.05 are marked with "X"
p_threshold <- 0.01
p_values <- ifelse(corr$p_adj<p_threshold, "X", "")

pheatmap(mat,
         breaks = breaks,
         color = colors,
         display_numbers = p_values,
         main = paste0("Correlations between bacteria and metabolites 
              (statistically significant associations (p < 0.05) marked with X)"),
         number_color = "yellow")


## -------------------------------------------------------------------------------------------------------------------------------
# Load package
library(biclust)

set.seed(19574)

# Find biclusters
bc <- biclust(corr$cor, method=BCPlaid(), fit.model = y ~ m, 
              background = TRUE, shuffle = 100, back.fit = 0, max.layers = 10, 
              iter.startup = 10, iter.layer = 100, verbose = FALSE)

bc


## -------------------------------------------------------------------------------------------------------------------------------
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


## ---- message=FALSE, warning=FALSE, fig.height=10, fig.width=16-----------------------------------------------------------------
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


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
reticulate::install_miniconda(force = TRUE)
reticulate::use_miniconda()
reticulate::py_install(packages = c("mofapy2"), pip = TRUE)


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
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


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
model_opts <- get_default_model_options(model)
model_opts$num_factors <- 15
head(model_opts)


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
train_opts <- get_default_training_options(model)
head(train_opts)


## ---- message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------
model.prepared <- prepare_mofa(
    object = model,
    model_options = model_opts
)
model.trained <- run_mofa(model.prepared)


## ---- message=FALSE, warning=FALSE, fig.height=8, fig.width=10------------------------------------------------------------------
library(patchwork)
wrap_plots(
    plot_variance_explained(model.trained, x="view", y="factor", plot_total = T),
    nrow = 2
) + plot_annotation(title = "Varience Explained per factor and assay",
                    theme = theme(plot.title = element_text(hjust = 0.5)))



## ---- warning=FALSE, message=FALSE, fig.height=8, fig.width=10------------------------------------------------------------------
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

butyrate_cutoff <- median(butyrate_df_train$Butyrate)
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
