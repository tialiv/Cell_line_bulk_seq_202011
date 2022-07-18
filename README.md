# Project 

Bulk RNA sequencing project on ovarian cell lines (COV434, KGN, primary ovarian cells) after exposure to DMSO, DES 10-6M & 10-10M and KTZ 10-5M & 10-9M.

## Analysis environment in R

Conda environment for analysis in R can be created by running the following code in the terminal using the environment.yml file

```

# Create conda environment
conda env create -f environment.yml -n rna

# Activate the environment
conda activate rna

# Running RStudio and run the code
rstudio &

# Deactivate the environment when finishing analysis
conda deactivate

# If package DOSE and clusterprofiler cannot be loaded and the problem is that DO.db needs to be reinstalled, activate the environment and run the following code
conda install -c bioconda bioconductor-do.db

```
