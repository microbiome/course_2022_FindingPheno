
# Getting started


## Checklist (before the workshop)

Install the following software in advance in order to avoid
unnecessary delays and leaving more time for the workshop contents.


* [R (version >4.1.0)](https://www.r-project.org/) 

* [RStudio](https://www.rstudio.com/products/rstudio/download/);
  choose "Rstudio Desktop" to download the latest version. Optional
  but preferred. For further details, check the [Rstudio home
  page](https://www.rstudio.com/).
  
* For Windows users: [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools40.html);
  Follow the instructions to install the toolkit. This might be required to compile some of the 
  packages required on this course.

* Install and load the required R packages


## Support and resources


For online support on installation and other matters, you can join us at:

 * Users: [miaverse Gitter channel](https://gitter.im/microbiome/miaverse?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
 * Developers: [Bioconductor Slack](https://bioc-community.herokuapp.com) #microbiomeexperiment channel (ask for an invitation)



## Installing and loading the required R packages

This section shows how to install and load all required packages into
the R session. Only uninstalled packages are installed.


```r
# List of packages that we need from cran and bioc 
packages <- c("BiocManager", "ggplot2", "pheatmap", "stringr", "igraph", "ANCOMBC",
             "microbiome", "httpuv", "microbiomeDataSets", "mia", "caret", "ranger",
            "dplyr", "miaViz", "knitr", "kableExtra", "vegan", "ecodist", "biclust",
            "patchwork", "pdp", "MLmetrics", "precrec")
```
 
The following script tries to load all required packages, and if they are not available, installs them.

```r
is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=T)
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)
```

```
##        BiocManager            ggplot2           pheatmap            stringr 
##               TRUE               TRUE               TRUE               TRUE 
##             igraph            ANCOMBC         microbiome             httpuv 
##               TRUE               TRUE               TRUE               TRUE 
## microbiomeDataSets                mia              caret             ranger 
##               TRUE               TRUE               TRUE               TRUE 
##              dplyr             miaViz              knitr         kableExtra 
##               TRUE               TRUE               TRUE               TRUE 
##              vegan            ecodist            biclust          patchwork 
##               TRUE               TRUE               TRUE               TRUE 
##                pdp          MLmetrics            precrec 
##               TRUE               TRUE               TRUE
```