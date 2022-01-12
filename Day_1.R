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