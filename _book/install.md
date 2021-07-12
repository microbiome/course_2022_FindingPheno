# Install and load required R packages

In this section, all required packages are installed and loaded into the
session. If packages are already installed, installation step is
skipped; only uninstalled packages are installed.

    # List of packages that we need from cran and bioc 
    cran_pkg <- c("BiocManager", "dplyr", "ecodist", "ggplot2", "gridExtra", "knitr", "vegan")
    bioc_pkg <- c("ANCOMBC", "ape", "DESeq2",  "DirichletMultinomial", "mia", "miaViz")

    # Gets those packages that are already installed
    cran_pkg_already_installed <- cran_pkg[ cran_pkg %in% installed.packages() ]
    bioc_pkg_already_installed <- bioc_pkg[ bioc_pkg %in% installed.packages() ]

    # Gets those packages that need to be installed
    cran_pkg_to_be_installed <- setdiff(cran_pkg, cran_pkg_already_installed)
    bioc_pkg_to_be_installed <- setdiff(bioc_pkg, bioc_pkg_already_installed)

    # If there are packages that need to be installed, installs them from CRAN
    if( length(cran_pkg_to_be_installed) ) {
       install.packages(cran_pkg_to_be_installed)
    }

    # If there are packages that need to be installed, installs them from Bioconductor
    if( length(bioc_pkg_to_be_installed) ) {
       BiocManager::install(bioc_pkg_to_be_installed, ask = F)
    }

Now all required packages are installed, so let’s load them into the
session. Some function names occur in multiple packages. That is why
miaverse’s packages mia and miaViz are prioritized. Packages that are
loaded first have higher priority.

    # Reorders bioc packages, so that mia and miaViz are first
    bioc_pkg <- c(bioc_pkg[ bioc_pkg %in% c("mia", "miaViz") ], bioc_pkg[ !bioc_pkg %in% c("mia", "miaViz") ] ) 

    # Loading all packages into session. Returns true if package was successfully loaded.
    sapply(c(bioc_pkg, cran_pkg), require, character.only = TRUE)

    ##                  mia               miaViz              ANCOMBC                  ape               DESeq2 DirichletMultinomial 
    ##                 TRUE                 TRUE                 TRUE                 TRUE                 TRUE                 TRUE 
    ##          BiocManager                dplyr              ecodist              ggplot2            gridExtra                knitr 
    ##                 TRUE                 TRUE                 TRUE                 TRUE                 TRUE                 TRUE 
    ##                vegan 
    ##                 TRUE