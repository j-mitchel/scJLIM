# scJLIM

Single-cell eQTL GWAS colocalization tool

## Installation

First install Julia:
```bash
# Linux
curl -fsSL https://install.julialang.org | sh
```

Then, restart your terminal:
```bash
exec $SHELL
```

Note: Do NOT install Julia via conda (conda install julia), as it may lead to incompatible system libraries and package loading errors.

Next, clone the scJLIM repo and create a conda environment with the necessary dependencies:
```bash
git clone https://github.com/j-mitchel/scJLIM.git
cd scJLIM
conda env create -f environment.yml
conda activate scjlim
# if using jupyter notebooks, also include the following
conda install -c conda-forge jupyter_client
```

Start R and install additional package dependencies:
```r
library(devtools)
install_deps(dependencies = TRUE)
install_github("yaowuliu/ACAT")
# if using jupyter notebooks, install IRkernel below
install.packages('IRkernel')
IRkernel::installspec(
name = "scjlim-r",
displayname = "R (scjlim)"
)
q()
```


## Walkthrough

The repository currently contains a [vignette](https://nbviewer.org/github/j-mitchel/scJLIM/blob/main/vignettes/tutorial.ipynb)
demonstrating how to use the package.

