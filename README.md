# Project 

Bulk RNA sequencing project on ovarian cell lines (COV434, KGN, primary ovarian cells) after exposure to DMSO, DES 10-6M & 10-10M and KTZ 10-5M & 10-9M.

## Analysis environment in R

Conda environment for analysis in R can be created using the environment.yml file

```

# Create conda environment
conda env create -f environment.yml -n rna

# Activate the environment
conda activate rna

# Running RStudio and run the code
rstudio &

# Deactivate the environment when finishing analysis
conda deactivate


project
|- doc/                documentation for the study
|
|- data/               raw and primary data, essentially all input files, never edit!
|  |- raw_external/
|  |- raw_internal/
|  |- meta/
|
|- code/               all code needed to go from input files to final results
|- notebooks/
|
|- intermediate/       output files from different analysis steps, can be deleted
|- scratch/            temporary files that can be safely deleted or lost
|- logs/               logs from the different analysis steps
|
|- results/            output from workflows and analyses
|  |- figures/
|  |- tables/
|  |- reports/
|
|- .gitignore          sets which parts of the repository that should be git tracked
|- Snakefile           project workflow, carries out analysis contained in code/
|- config.yml          configuration of the project workflow
|- environment.yml     software dependencies list, used to create a project environment
|- Dockerfile          recipe to create a project container
```
