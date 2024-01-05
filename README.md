# SNACS

Single Nucleotide Polymorphism (SNP) and Antibody-Based Cell Sorting (SNACS) is an R package for demultiplexing single-cell DNA sequencing data using DNA sequencing data and hash antibody data.

'SNACS' package can be installed as follows:

Launch R.

'SNACS' is dependent on 'heatmap4'. 'heatmap4' is available in r-universe. In order to automatically install 'heatmap4' with 'SNACS', the repository of the former package has to be specified before installing 'SNACS'.
> options(repos = c(
    rituroy = 'https://rituroy.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

Then install 'SNACS'.
> remotes::install_github('olshena/SNACS',subdir='/SNACS',build_vignettes=TRUE)

'heatmap4' can also be installed directly from 'rituroy' universe.
> install.packages('heatmap4', repos = 'https://github.com/rituroy/heatmap4')

Load 'SNACS'.
> require(SNACS)

See a demo of 'SNACS'.
> vignette('SNACS')
