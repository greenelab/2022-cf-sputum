#devtools::install_github("greenelab/TDM")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager") # should be installed in conda environment already

BiocManager::install("greenelab/ADAGEpath", update = FALSE) # will install dependency greenelab/TDB automatically

file.create(snakemake@output[['install']])
