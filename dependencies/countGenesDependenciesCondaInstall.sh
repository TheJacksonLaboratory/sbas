#!/bin/bash
conda install mamba -y
mamba install -c bioconda bioconductor-genomeinfodb -y
mamba install -c conda-forge r-vctrs -y
mamba install -c conda-forge r-rlang -y
mamba install -c r r-dplyr -y
mamba install -c bioconda bioconductor-biobase -y
mamba install -c conda-forge r-tibble -y
mamba install -c conda-forge r-r.utils -y
mamba install -c bioconda bioconductor-rtracklayer -y
