# scJLIM

Single-cell eQTL GWAS colocalization tool

## Installation
To install the latest version of scITD from GitHub:

``` r
# Install the required JuliaCall package and Julia language for mixed model fitting
install.packages("JuliaCall")
library(JuliaCall)
install_julia()

# Install the required ACAT package for calculating cauchy combined p-values
library(devtools)
install_github("yaowuliu/ACAT")

# Install scJLIM package
install_github("j-mitchel/scJLIM")
```

## Walkthrough
The repository currently contains a [vignette](https://github.com/j-mitchel/scJLIM/blob/main/vignettes/tutorial.ipynb)
demonstrating how to use the package.

