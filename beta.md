# Beta diversity

This notebook shows how to analyse and visualize beta diversity.

Beta diversity reflects the difference in microbial composition between
two samples. Similar samples have a low beta diversity.

Several different distance metrics are available. Some of the common
choices include Bray-Curtis, Unifrac, Jaccard, and Aitchison index. Each
of these dissimilarity measures emphasize different aspects of
similarity.

## Examples of PCoA with different settings

After estimating beta diversity we can vizualize sample similarity with
dimension reduction techniques such as Principal Coordinate Analysis
(PCoA).

PCoA takes a dissimilarity matrix as input. The output is usually a 2 or
3-dimensional euclidean space. The idea is to project the data so that
the distances between different samples are maximized. The projection is
often non-linear and designed to reveal local or global structures in
the data distribution.

### PCoA for ASV-level data with Bray-Curtis

Let us next show how to visualize sample similarities using the PCoA
method and a selected dissimilarity measure. We use the same tse data
object defined in the earlier notebooks.

    # Relative abundance table
    rel_abund_assay <- assays(tse)$relabundance

    # Transposes it to get taxa to columns
    rel_abund_assay <- t(rel_abund_assay)

    # Calculates Bray-Curtis distances between samples. Because taxa is in columns,
    # it is used to compare different samples.
    bray_curtis_dist <- vegan::vegdist(rel_abund_assay, method = "bray")

    # Does principal coordinate analysis
    bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)

    # Creates a data frame from principal coordinates
    bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                      pcoa2 = bray_curtis_pcoa$vectors[,2])

    # Creates a plot
    bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
      geom_point() +
      labs(x = "Coordinate 1",
           y = "Coordinate 2", 
           title = "Bray-Curtis PCoA with relative abundances") 
      theme(title = element_text(size = 10)) # makes titles smaller

    ## List of 1
    ##  $ title:List of 11
    ##   ..$ family       : NULL
    ##   ..$ face         : NULL
    ##   ..$ colour       : NULL
    ##   ..$ size         : num 10
    ##   ..$ hjust        : NULL
    ##   ..$ vjust        : NULL
    ##   ..$ angle        : NULL
    ##   ..$ lineheight   : NULL
    ##   ..$ margin       : NULL
    ##   ..$ debug        : NULL
    ##   ..$ inherit.blank: logi FALSE
    ##   ..- attr(*, "class")= chr [1:2] "element_text" "element"
    ##  - attr(*, "class")= chr [1:2] "theme" "gg"
    ##  - attr(*, "complete")= logi FALSE
    ##  - attr(*, "validate")= logi TRUE

    bray_curtis_plot

![](beta_files/figure-markdown_strict/pcoa_asv_bc-1.png)

### PCoA for ASV-level data with Aitchison distance

Aitchison distance corresponds to Euclidean distances between CLR
transformed sample abundance vectors.

    # Does clr transformation. Pseudocount is added, because data contains zeros. 
    tse <- transformCounts(tse, method = "clr", pseudocount = 1)

    # Gets clr table
    clr_assay <- assays(tse)$clr

    # Transposes it to get taxa to columns
    clr_assay <- t(clr_assay)

    # Calculates Euclidean distances between samples. Because taxa is in columns,
    # it is used to compare different samples.
    euclidean_dist <- vegan::vegdist(clr_assay, method = "euclidean")


    # Does principal coordinate analysis
    euclidean_pcoa <- ecodist::pco(euclidean_dist)

    # Creates a data frame from principal coordinates
    euclidean_pcoa_df <- data.frame(pcoa1 = euclidean_pcoa$vectors[,1], 
                                    pcoa2 = euclidean_pcoa$vectors[,2])

    # Creates a plot
    euclidean_plot <- ggplot(data = euclidean_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
      geom_point() +
      labs(x = "Coordinate 1",
           y = "Coordinate 2",
           title = "Euclidean PCoA with CLR transformation") +
      theme(title = element_text(size = 12)) # makes titles smaller

    euclidean_plot

![](beta_files/figure-markdown_strict/pcoa_asv_aitchison-1.png)

### PCoA aggregated to Phylum level

We use again the Aitchison distances in this example but this time we
use data that was aggregated to the phylum level in the earlier
examples.

    # Does clr transformation. Psuedocount is added, because data contains zeros. 
    tse_phylum <- transformCounts(tse_phylum, method = "clr", pseudocount = 1)

    # Gets clr table
    clr_phylum_assay <- assays(tse_phylum)$clr

    # Transposes it to get taxa to columns
    clr_phylum_assay <- t(clr_phylum_assay)

    # Calculates Euclidean distances between samples. Because taxa is in columns,
    # it is used to compare different samples.
    euclidean_phylum_dist <- vegan::vegdist(clr_assay, method = "euclidean")

    # Does principal coordinate analysis
    euclidean_phylum_pcoa <- ecodist::pco(euclidean_phylum_dist)

    # Creates a data frame from principal coordinates
    euclidean_phylum_pcoa_df <- data.frame(pcoa1 = euclidean_phylum_pcoa$vectors[,1], 
                                           pcoa2 = euclidean_phylum_pcoa$vectors[,2])

    # Creates a plot
    euclidean_phylum_plot <- ggplot(data = euclidean_phylum_pcoa_df, aes(x=pcoa1, y=pcoa2)) +
      geom_point() +
      labs(x = "Coordinate 1",
           y = "Coordinate 2",
           title = "Aitchison distances at Phylum level") +  
      theme(title = element_text(size = 12)) # makes titles smaller

    euclidean_phylum_plot

![](beta_files/figure-markdown_strict/pcoa_phylum_aitchison-1.png)

## Highlighting external variables on PCoA plot

### PCoA with discrete sample grouping variable shown with colors

We can add grouping variable to existing plots. Let’s add coloring to
the CLR transformed, Genus level PCoA.

    # Adds coloring information to the data frame, creates new column
    euclidean_patient_status_pcoa_df <- cbind(euclidean_pcoa_df,
                                 patient_status = colData(tse)$patient_status)

    # Creates a plot
    euclidean_patient_status_plot <- ggplot(data = euclidean_patient_status_pcoa_df, 
                                            aes(x=pcoa1, y=pcoa2,
                                                color = patient_status)) +
      geom_point() +
      labs(x = "Coordinate 1",
           y = "Coordinate 2",
           title = "PCoA with Aitchison distances") +
      theme(title = element_text(size = 12)) # makes titles smaller

    euclidean_patient_status_plot

![](beta_files/figure-markdown_strict/pcoa_genus-1.png)

### PCoA with continuous sample grouping variable shown with colors

We can also use continues values to group variables. Let’s use Shannon
diversity index that we calculated in “Alpha diversity” notebook.

    # Adds coloring information to the data frame, creates new column
    euclidean_shannon_pcoa_df <- cbind(euclidean_pcoa_df,
                                 shannon = colData(tse)$Shannon_index)

    # Creates a plot
    euclidean_shannon_plot <- ggplot(data = euclidean_shannon_pcoa_df, 
                                     aes(x=pcoa1, y=pcoa2,
                                         color = shannon)) + 
      geom_point() +
      labs(x = "Coordinate 1",
           y = "Coordinate 2",
           title = "PCoA with Aitchison distances") +
      theme(title = element_text(size = 12)) # makes titles smaller

    euclidean_shannon_plot

![](beta_files/figure-markdown_strict/pcoa_coloring-1.png)

## Estimating associations with an external variable

Next we show how the quantify the strength of association between the
variation in community composition (beta diversity) and external
factors.

The standard way to do this is to perform a so-called permutational
multivariate analysis of variance (PERMANOVA) test.

    # Relative abundance table
    rel_abund_assay <- assays(tse)$relabundance

    # Transposes it to get taxa to columns
    rel_abund_assay <- t(rel_abund_assay)

    permanova_cohort <- vegan::adonis(rel_abund_assay ~ cohort,
                                      data = colData(tse),
                                      permutations = 9999)

    # P-value
    print(paste0("Different different cohorts and variance of abundance between samples, p-value: ", 
                 as.data.frame(permanova_cohort$aov.tab)["cohort", "Pr(>F)"]))

    ## [1] "Different different cohorts and variance of abundance between samples, p-value: 0.7477"

As we see, the cohort variable is not significantly associated with
microbiota composition (p-value is over 0.05).

We can, however, visualize those taxa whose abundance are the most
different between cohorts. This gives us information which taxa’s
abundances tend to differ between different cohorts.

In order to do that, we need coefficients of taxa.

    # Gets the coefficients
    coef <- coefficients(permanova_cohort)["cohort1",]

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

![](beta_files/figure-markdown_strict/unnamed-chunk-2-1.png)

The above plot shows taxa as code names, and it is hard to tell which
bacterial groups they represent. However, it is easy to add human
readable names. We can fetch those from rowData. Here we use Genus level
names.

    # Gets corresponding Genus level names and stores them to top.coef
    names <- rowData(tse)[names(top.coef), ][,"Genus"]

    # Adds new labels to the plot
    top_taxa_coeffient_plot <- top_taxa_coeffient_plot +
      scale_y_discrete(labels = names) # Adds new labels

    top_taxa_coeffient_plot

![](beta_files/figure-markdown_strict/unnamed-chunk-3-1.png)

## Further resources

For more examples, see a dedicated section on beta diversity in the
[online book](https://microbiome.github.io/OMA/).