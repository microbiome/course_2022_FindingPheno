# Introduction to multi-omics data analysis 

**Welcome to the multi-omics data analysis workshop**

<img src="https://user-images.githubusercontent.com/60338854/121848694-1072a480-ccf3-11eb-9af2-7fdefd8d1794.png" alt="ML4microbiome" width="50%"/>

Figure source: Moreno-Indias _et al_. (2021) [Statistical and Machine Learning Techniques in Human Microbiome Studies: Contemporary Challenges and Solutions](https://doi.org/10.3389/fmicb.2021.635781). Frontiers in Microbiology 12:11. 


## Rendering the book

You can render the book locally in R with:

```{r serve}
bookdown::serve_book()
``` 

## The miaverse framework

The [_miaverse_](https://microbiome.github.io) (mia = **MI**crobiome **A**nalysis) is an
R/Bioconductor framework for microbiome data science. It aims to
extend the capabilities of another popular framework,
[phyloseq](https://joey711.github.io/phyloseq/).

The miaverse framework consists of an efficient data structure, an
associated package ecosystem, demonstration data sets, and open
documentation. These are explained in more detail in the online book
[Orchestrating Microbiome Analysis](https://microbiome.github.io/OMA).

This workshop material walks you through example workflows for multi-omics data
analysis covering data access, exploration, analysis, visualization and reproducible
reporting. **You can run the workflow by simply copy-pasting the
examples.** For advanced material, you can test and modify further
examples from the [OMA book](https://microbiome.github.io/OMA), or try
to apply the techniques to your own data.




# Learning goals [TO DO]

This workshop provides an overview of bioinformatics
tools for multi-omics studies, ranging from data
preprocessing to statistical analysis and reproducible reporting.


**Target audience** Advanced students and applied researchers who wish
  to develop their skills in microbial community analysis. [TO DO]

**Venue** [TO DO]





## Acknowledgments

**Citation** "Introduction to miaverse (2021). Tuomas Borman, Henrik Eckermann, Chouaib Benchraka, Leo Lahti. URL: https://microbiome.github.io".

**Contact**
- [Leo Lahti](http://datascience.utu.fi), University of Turku
- [mia Collective](https://microbiome.github.io)

**License** All material is released under the open [CC BY-NC-SA 3.0 License](LICENSE).

- Landing page (html): [workshop teaching material](https://microbiome.github.io/course_2022_FindingPheno/)
- Source code (github): [workshop teaching material](https://github.com/microbiome/course_2022_FindingPheno)

The source code of this repository is fully reproducible and contains
the Rmd files with executable code. All files can be rendered at one
go by running the file [main.R](main.R). You can check the file for
details on how to clone the repository and convert it into a gitbook,
although this is not necessary for the training.
