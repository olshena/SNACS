# SNACS

Single Nucleotide Polymorphism (SNP) and Antibody-Based Cell Sorting (SNACS) is an R package for demultiplexing single-cell DNA sequencing data using DNA sequencing data and hash antibody data.

# Installation

Launch R.

'SNACS' is dependent on 'heatmap4'. 'heatmap4' is available in r-universe. In order to automatically install 'heatmap4' with 'SNACS', the repository of the former package has to be specified before installing 'SNACS'.
```{r}
options(repos = c(
    rituroy = 'https://rituroy.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))
```

Then install 'SNACS'.
```{r}
remotes::install_github('olshena/SNACS',subdir='/SNACS',build_vignettes=TRUE)
```

If you would like to run SNACS using an hdf5 file as an input, you will also need the hdf5 package from Bioconductor. 
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
```

# Usage
Load 'SNACS'.
```{r}
library('SNACS')
```

See a demo of 'SNACS'.
```{r}
vignette('SNACS')
```
