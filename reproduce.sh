# This script can be only executed on machines with conda already available.
# Please make sure you have initialised the terminal for use of conda before proceeding
# To do so you can run the following command in your terminal:
# `conda init zsh && exec -l zsh`
# or
#  `conda init bash && exec -l bash`


# Clone sbas repo
git clone https://github.com/TheJacksonLaboratory/sbas

# cd into repo
cd sbas

# Install dependencies in your linux machine with conda available
## Install mamba, a faster alternative/implementation compared conda 
conda install mamba -y

## Create a new isolated environment for the analysis
mamba env create --name sbas -f environment.yml 

## Activate the new environment
conda activate sbas

# Retrieve prerequisite input files for Jupyter Notebooks from ZENODO
wget https://zenodo.org/record/4179559/files/as.tar.gz
wget https://zenodo.org/record/4179559/files/dge.tar.gz
wget https://zenodo.org/record/4179559/files/fromGTF.tar.gz
wget https://zenodo.org/record/4179559/files/gtex.tar.gz 
wget https://zenodo.org/record/4179559/files/rmats_final.tar.gz
wget https://zenodo.org/record/4179559/files/srr.tar.gz

# Decompress archives into the empty data folder and delete the archives after 
tar xzvf as.tar.gz -C data && rm as.tar.gz
tar xzvf dge.tar.gz -C data && rm dge.tar.gz
tar xzvf fromGTF.tar.gz -C data && rm fromGTF.tar.gz
tar xzvf gtex.tar.gz  -C data && rm  gtex.tar.gz
tar xzvf rmats_final.tar.gz -C data && rm rmats_final.tar.gz
tar xzvf srr.tar.gz -C data && rm srr.tar.gz

# cd into jupyter
cd jupyter

# Execute programmatically the notebooks with Papermill
papermill countGenesAndEvents.ipynb countGenesAndEvents.ipynb
papermill expressionHeatplot.ipynb expressionHeatplot.ipynb 
papermill totalDGEByTissue.ipynb totalDGEByTissue.ipynb 
papermill alternativeSplicingHeatplot.ipynb alternativeSplicingHeatplot.ipynb 
papermill totalAlternativeSplicingByTissue.ipynb totalAlternativeSplicingByTissue.ipynb
papermill XchromosomalEscape.ipynb XchromosomalEscape.ipynb
papermill splicingIndex.ipynb splicingIndex.ipynb
papermill spliceTypeByChromosome.ipynb spliceTypeByChromosome.ipynb
papermill altSplicing_events_per_gene.ipynb altSplicing_events_per_gene.ipynb
papermill tissue_piechart.ipynb tissue_piechart.ipynb
