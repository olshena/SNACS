# SNACS

SNACS is a package in R. It demuliplexes samples of single cell mutation data using single cell mutation data and hash antibody data.

To install "SNACS" package we will need "devtools" package.
"SNACS" is dependent on "heatmap4" which is available on GitHub.

> library(devtools)

Run this to install "heatmap4".
> install_github("UCSF-CBI/heatmap4/heatmap4")
or 
> install_github("UCSF-CBI/heatmap4/heatmap4",auth_token="abc")

> install_github("olshena/SNACS/SNACS",build_vignettes=TRUE)
or
> install_github("olshena/SNACS/SNACS",build_vignettes=TRUE,auth_token="abc")

Run this to install "SNACS".
> install_github("olshena/SNACS/SNACS",build_vignettes=TRUE)
or
> install_github("olshena/SNACS/SNACS",build_vignettes=TRUE,auth_token="abc")

This will give a demo of the package.
> require(SNACS)
or
> vignette("SNACS")
