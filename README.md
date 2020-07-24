# The impact of sex on alternative splicing

This repository documents the analysis performed for [The impact of sex on alternative splicing](https://www.biorxiv.org/content/10.1101/490904v1.full);
note that a manuscript with a modified version of the analysis has been submitted. To reproduce the analysis, users will need to go through several steps.

#. Get access to the [Genotype-Tissue Expression (GTEx)](https://www.gtexportal.org/home/) RNAseq data (an application to dbGAP for access to the dataset phs000424.v8.v2 is required)
#. Align each RNAseq sample using hisat2 and create a matrix of counts for each of a variety of splicing types was generated by the rMATS. See HERE for a nextflow script. The script can be run on the [cloudOS/lifebit platform](https://lifebit.ai/).
#. Run the Jupyter notebooks from this repository to perform the individual analyses.


This repository documents the interactive analysis for the results of running the rmats-nf pipeline.

## 1. Get access

The RNA-seq samples analyzed in this project are restricted access (dbGAP phs000424.v8.v2). See the
[database of Genotypes and Phenotypes (dbGaP)](https://www.ncbi.nlm.nih.gov/gap/) for details.

## 2. Processing the RNA-seq samples

See the manuscript for methods details. In brief, we ran the nextflow script at https://github.com/lifebit-ai/rmats-nf  to align the RNA-seq samples with hisat2 and to characterize splicing events with rMATS. Results from individual samples are summarized in 'matrix' files. To run the Jupyter scripts in the
next section, you will need to place these files in a results bucket (if you are using the cloudos system) or in some other defined location.

## 3. Running the notebooks

Each of the results described in the manuscript was generated by one or more Jupyter notebooks in this repository.
There are a number of R packages that need to be installed prior to running the notebooks. This process is described from the
cloudos environment in [this document](https://github.com/TheJacksonLaboratory/sbas/blob/master/SettingUpRobinsonLabNotebook.MD). If running the notebooks in another environment, simply run the 
setup scripts. 

### 3.1 Summarizing events

Most of the notebooks require that the raw rMATS files are first processed to generate summary files. This is done by the notebook
[countGenesAndEvents.ipynb](https://github.com/TheJacksonLaboratory/sbas/blob/master/countGenesAndEvents.ipynb). The following list summarizes
what each of the notebooks does.

#. [countGenesAndEvents.ipynb](https://github.com/TheJacksonLaboratory/sbas/blob/master/countGenesAndEvents.ipynb). Set up the overall analysis. Write various files to the ``data`` subdirectory that will be used by other scripts.
#. [XchromosomalEscapePlot.ipynb](https://github.com/TheJacksonLaboratory/sbas/blob/master/XchromosomalEscapePlot.ipynb). Investigate the overlap of alternative splicing and genes on the X chromosome that escape inactivation.


