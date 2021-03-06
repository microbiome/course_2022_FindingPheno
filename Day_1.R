################################################################################
#####**********************************************************************#####
#####*********************** DAY 1 PRACTICAL ******************************#####
#####**********************************************************************#####
################################################################################

################################################################################
# Chapter 4 Data
################################################################################

## 4.3 Importing data in R

library(stringr)

# Load the data
mae <- microbiomeDataSets::HintikkaXOData()

# Drop off those bacteria that do not include information in Phylum or lower levels
mae[[1]] <- mae[[1]][!is.na(rowData(mae[[1]])$Phylum), ]

# Clean taxonomy data, so that names do not include addtional characters
rowData(mae[[1]]) <- DataFrame(apply(rowData(mae[[1]]), 2,
                                     str_remove, pattern = "._[0-9]__"))

mae

################################################################################
# Chapter 5 Microbiome data exploration
################################################################################

## 5.1 Data structure

# mae[[1]]: indexing/retrieving the taxonomic data experiment
dim(mae[[1]])
head(rowData(mae[[1]]))

# For simplicity, classify all high-fat diets as high-fat, and all the low-fat
# diets as low-fat diets
colData(mae)$Diet <- ifelse(colData(mae)$Diet == "High-fat" |
                                colData(mae)$Diet == "High-fat + XOS",
                            "High-fat", "Low-fat")

head(colData(mae))

sort(table(colData(mae)$Diet))

### 5.1.1 Transformations

# Calculates relative abundances, and stores the table to assays
mae[[1]] <- transformCounts(mae[[1]], method = "relabundance")

### 5.1.2 Aggregation

se_phylum <- agglomerateByRank(mae[[1]], rank = "Phylum")

# Show dimensionality
dim(se_phylum)

head(rowData(se_phylum))

temp <- rowData(agglomerateByRank(mae[[1]], rank = "Genus"))

# Prints those taxa that do not have information at the Genus level (NA)
head(temp[which(is.na(temp$Genus)),])

temp2 <- rowData(agglomerateByRank(mae[[1]], rank = "Genus", na.rm = TRUE))
print(paste0("Agglomeration with na.rm = FALSE: ", dim(temp)[1], " taxa."))
print(paste0("Agglomeration with na.rm = TRUE: ", dim(temp2)[1], " taxa."))

## 5.2 Visualization

# Here we specify "relabundance" to be abundance table that we use for plotting.
# Note that we can use agglomerated or non-agglomerated mae[[1]] as an input, because
# the function agglomeration is built-in option.

# Legend does not fit into picture, so its height is reduced.
plot_abundance <- plotAbundance(mae[[1]], abund_values="relabundance", rank = "Phylum") +
    theme(legend.key.height = unit(0.5, "cm")) +
    scale_y_continuous(label = scales::percent)

plot_abundance

# Subset data by taking only Firmicutes
se_firmicutes <- se_phylum["Firmicutes"]

# Gets the abundance table
abundance_firmicutes <- assay(se_firmicutes, "relabundance")

# Creates a data frame object, where first column includes abundances
firmicutes_abund_df <- as.data.frame(t(abundance_firmicutes))
# Rename the first and only column
colnames(firmicutes_abund_df) <- "abund"

# Creates a plot. Parameters inside feom_density are optional. With
# geom_density(bw=1000), it is possible to adjust bandwidth.
firmicutes_abund_plot <- ggplot(firmicutes_abund_df, aes(x = abund)) +
    geom_density(color="darkred", fill="lightblue") +
    labs(x = "Relative abundance", title = "Firmicutes") +
    theme_classic() + # Changes the background
    scale_x_continuous(label = scales::percent)

firmicutes_abund_plot


################################################################################
# Chapter 8 Beta diversity
################################################################################

# Beta diversity is another name for sample dissimilarity. It quantifies
# differences in the overall taxonomic composition between two samples.
#
# Common indices include Bray-Curtis, Unifrac, Jaccard index, and the
# Aitchison distance. For more background information
# and examples, you can check the dedicated section in online book


# ## Examples of PCoA with different settings
#
# Beta diversity estimation generates a (dis)similarity matrix that
# contains for each sample (rows) the dissimilarity to any other sample
# (columns).
#
# This complex set of pairwise relations can be visualized in
# informative ways, and even coupled with other explanatory
# variables. As a first step, we compress the information to a lower
# dimensionality, or fewer principal components, and then visualize
# sample similarity based on that using ordination techniques, such as
# Principal Coordinate Analysis (PCoA). PCoA is a non-linear dimension
# reduction technique, and with Euclidean distances it is is identical
# to the linear PCA (except for potential scaling).
#
# We typically retain just the two (or three) most informative top
# components, and ignore the other information. Each sample has a score
# on each of these components, and each component measures the variation
# across a set of correlated taxa. The top components are then easily
# visualized on a two (or three) dimensional display.
#
# Let us next look at some concrete examples.
#
# ### PCoA for ASV-level data with Bray-Curtis
#
# Let us start with PCoA based on a Bray-Curtis dissimilarity matrix.

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


# ## Highlighting external variables
#
# We can map other variables on the same plot for example by coloring
# the points accordingly.
#
# The following is an example with a discrete grouping variable (Diet) shown with colors:
#

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
#
# ## Estimating associations with an external variable
#
# Next to visualizing whether any variable is associated with
# differences between samples, we can also quantify the strength of the
# association between community composition (beta diversity) and
# external factors.
#
# The standard way to do this is to perform a so-called permutational
# multivariate analysis of variance (PERMANOVA). This method takes as
# input the abundance table, which measure of distance you want to base
# the test on and a formula that tells the model how you think the
# variables are associated with each other.


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
#
#
# The diet variable is significantly associated with
# microbiota composition (p-value is less than 0.05).
#
# We can visualize those taxa whose abundances drive the
# differences between diets. We first need to extract the model
# coefficients of taxa:


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

#
# The above plot shows taxa as code names, and it is hard to tell which
# bacterial groups they represent. However, it is easy to add human readable
# names. We can fetch those from our rowData. Here we use Genus level names:


# Gets corresponding Genus level names and stores them to top.coef
names <- rowData(mae[[1]])[names(top.coef), ][,"Genus"]

# Adds new labels to the plot
top_taxa_coeffient_plot <- top_taxa_coeffient_plot +
    scale_y_discrete(labels = names) # Adds new labels
top_taxa_coeffient_plot


# The same test can be conducted using the ordination from PCoA as follows:


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
    labs(title = paste0("PCoA beta diversity ordination for microbiome samples"), x = paste0("PC1 (p = ", p_values[["pcoa1"]], ")"), y = paste0("PC2 (p = ", p_values[["pcoa2"]], ")")) +
    theme_bw()

plot

#
# There are many alternative and complementary methods for analysing
# community composition. For more examples, see a dedicated section on
# # beta diversity in the online book
#
#  ## Community typing
#
#
#  A dedicated section presenting examples on community typing is in the
#  online book
