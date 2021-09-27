# This script can be only executed on machines with conda already available.
# Please make sure you have initialised the terminal for use of conda before proceeding
# To do so you can run the following command in your terminal:
# `conda init zsh && exec -l zsh`
# or
#  `conda init bash && exec -l bash`
# Clone sbas repo
# git clone https://github.com/TheJacksonLaboratory/sbas
# cd into repo
# cd sbas

# Install dependencies in your linux machine with conda available

## Create new empty environment to avoid conflicts
# prior to running this script please run
#
#`conda create --name sbas -y `
#
# and then Activate the new environment
#
# with
# `conda activate sbas`
#
## Execute this script now with
#
# `bash reproduce.sh`
#
## Install mamba, a faster alternative/implementation compared conda
conda install mamba -y

# install github command line interface for convenience
mamba install gh --channel conda-forge -y

# install gcsfs due to a bug noted in issue on github 
mamba install  gcsfs==0.2.3 --force-reinstall -y

# install papermill inside our isolated environment
mamba install papermill -y

# update the isolated environment for the analysis
mamba env update --name sbas -f environment.yml

# install the interactive r kernel for the jupyterlab notebooks
mamba install -c r r-irkernel -y

# install the interactive pythong kernel for the jupyterlab notebooks
mamba install -c conda-forge ipykernel -y

# intall the specific notebook dependencies
source dependencies/alternativeSplicingHeatmapDependenciesCondaInstall.sh
source dependencies/countGenesDependenciesCondaInstall.sh
source dependencies/differentialGeneAnalysisDependenciesCondaInstall.sh
source dependencies/differentialSplicingJunctionAnalysis.sh
source dependencies/spliceTypeByChromosomeDependenciesCondaInstall.sh
source dependencies/splicingIndexDependenciesCondaInstall.sh

# Retrieve prerequisite input files for Jupyter Notebooks from ZENODO
wget https://zenodo.org/record/5524975/files/fromGTF.tar.gz
wget https://zenodo.org/record/5524975/files/gtex.tar.gz 
wget https://zenodo.org/record/5524975/files/rmats_final.tar.gz
wget https://zenodo.org/record/5524975/files/srr.tar.gz
wget https://zenodo.org/record/5524975/files/SraRunTable.txt.gz
# Decompress archives into the empty data folder and delete the archives after 
tar xzvf fromGTF.tar.gz -C data && rm fromGTF.tar.gz
tar xzvf gtex.tar.gz  -C data && rm  gtex.tar.gz
tar xzvf rmats_final.tar.gz -C data && rm rmats_final.tar.gz
tar xzvf srr.tar.gz -C data && rm srr.tar.gz
# we do not unzip SraRunTable.txt.gz - it is used in its gz state

# cd into jupyter
cd jupyter

# Execute programmatically the notebooks with Papermill
papermill differentialGeneExpressionAnalysis.ipynb differentialGeneExpressionAnalysis.ipynb -k ir
papermill differentialSplicingJunctionAnalysis.ipynb differentialSplicingJunctionAnalysis.ipynb -k ir
papermill countGenesAndEvents.ipynb countGenesAndEvents.ipynb -k ir
papermill expressionHeatplot.ipynb expressionHeatplot.ipynb -k ir
papermill totalDGEByTissue.ipynb totalDGEByTissue.ipynb -k ir
papermill alternativeSplicingHeatplot.ipynb alternativeSplicingHeatplot.ipynb -k ir
papermill totalAlternativeSplicingByTissue.ipynb totalAlternativeSplicingByTissue.ipynb -k ir
papermill XchromosomalEscape.ipynb XchromosomalEscape.ipynb -k ir
papermill splicingIndex.ipynb splicingIndex.ipynb -k ir
papermill spliceTypeByChromosome.ipynb spliceTypeByChromosome.ipynb -k ir
papermill altSplicing_events_per_gene.ipynb altSplicing_events_per_gene.ipynb -k ir
papermill tissue_piechart.ipynb tissue_piechart.ipynb -k ir
