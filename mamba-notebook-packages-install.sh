#!/bin/bash
# environment.yml does not seem to update the r-packages at this time
#
mamba install -c conda-forge ipykernel -y
mamba install -c r r-irkernel -y
mamba install -c bioconda bioconductor-biobase -y
mamba install -c bioconda bioconductor-rtracklayer -y
mamba install -c bioconda bioconductor-edger -y
mamba install -c conda-forge nb_conda_kernels -y
mamba install -c r r-icestaf -y
