# Getting set up fast on the Robinson-Lab Jupyter Lab notebook.

## 1 Cloning the repository

You should begin by git cloning the repository
At the launcher page

### 1.1 Start a terminal

From the Launcher Page select Terminal

![Launcher Page](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/LauncherPage.png)

### 1.1 Clone the repository

Within the terminal window, git clone the Robinson Lab sbas repository

```bash
git clone https://github.com/TheJacksonLaboratory/sbas.git
```

![git clone sbas](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/gitCloneRepository.png)

### 1.2. Get your files

The specific hash will be different for each mount.  The filename completion by using tab will allow you to quickly get your files.

Only `two` of the mounted files are really needed to get your work done.

#### 1.2.1 data.tar.gz

This command will unpack the data.tar.gz into your sbas/data directory

```bash
cd sbas
tar xvzf ../mounted_data/{long hashed key}data.tar.gz{more hashed data} -C data
```

Which looks like this in my window on Peters notebook.

![tar xvzf data.tar.gz](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/unpackingDataTarAndContents.png)

#### 1.2.2 DGE_splicing_data.tar.gz

```bash
tar xvzf ../mounted_data/{long hashed key}DGE_splicing_data.tar.gz{more hashed data} -C data
```

![tar xvzf DGE_splicing_data.tar.gz](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/unpackingFilesFromMountedDirectory.png)

## 2 Preparing your jupyter notebook

After getting your notebook started, you need to prepare for your notebooks expectations regarding librarys it is interested in using.  To do so, do the following in a terminal window.
Conda install is your friend.   Saying yes at the end means you can fire and forget.   Ones says yes with "-y" at the end of the conda install command.

Before opening your notebook - set up the minimal environment -- you may need more than this but this will get you started to run `countGenesAndEvents.ipynb` for example.

### 2.1 These conda install commands below get you going to start. 

```bash
conda install -c bioconda bioconductor-edger -y
conda install -c conda-forge r-readr -y
conda install -c bioconda bioconductor-multtest -y
conda install -c conda-forge r-r.utils -y 
conda install -c conda-forge r-snakecase -y
conda install -c conda-forge r-statmod -y
conda install -c bioconda bioconductor-stan -y
conda install -c bioconda bioconductor-rtracklayer -y
conda install -c conda-forge r-reshape -y
conda install -c conda-forge r-pheatmap -y
conda install -c bioconda bioconductor-limma -y
conda install -c bioconda bioconductor-edger -y
conda install -c bioconda bioconductor-multtest -y
conda install -c conda-forge r-r.utils -y
conda install -c bioconda bioconductor-rtracklayer -y
conda install -c conda-forge r-ggsci -y
conda install -c conda-forge r-rjags -y
conda install -c conda-forge r-compute.es -y
conda install -c conda-forge plotly -y
conda install -c bioconda r-pheatmap -y
conda install -c bioconda gffread -y
```

### 2.2 What do you do if there are more missing packages

Refer to the documentation to go further with Christina’s gitbook.

https://lifebit.gitbook.io/cloudos/web-interface/jupyter-lab#how-can-i-install-the-dependencies-i-need-in-the-notebook-session

## 3 Opening your notebook and unpacking the data files.

It is worthwhile that everyone run this notebook `countGenesAndEvents.ipynb` first as it will unpack the tarred files and show you the code that can read them in.  For any figure work you can use this as a template to create your work.

Now to test that all is well start with this notebook.

`countGenesAndEvents.ipynb`

Running this notebook will give some of the required downstream files for the subsequent notebooks.

### 3.1 Click Folder Icon

Clicking the folder icon in the upper left of the workspace will show you the mounted and clone repository directories.

![click folder icon](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/ClickFolderIconWillShowYoutheMountedAndClonedRepository.png)

### 3.2 Click sbas folder to start

Clicking the `sbas` folder will show you the subdirectories, including the location of the jupyter notebooks.

![click sbas folder](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/ClickingSBASWillShowSubDirectoriesAndJupyter.png)

### 3.3 Click jupyter folder

Clicking the `jupyter` folder will show you the jupyter notebooks we have.

![click sbas folder](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/ClickingJupyterWIllShowTheNotebooks.png)

### 3.4 Click countGenesAndEvents.ipynb to start

Clicking the first notebook to run, `countGenesAndEvents.ipynb` will show you the unpacking of data files and how to upload them into your jupyter notebook to begin your analyses.

![click sbas folder](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/ClickCountGenesAndEventsToStartAndUnloadAndUploadVariables.png)

### 3.5 Restart Kernal and Clear All Outputs

To have a clean start within this environment, click `kernal` and select `Interupt Kernal and clear all outputs` to have a clean and fresh start.

![click sbas folder](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/AGoodThingToDoToHaveACleanSlateIsRestartAndClearOutputs.png)

### 3.6 Run the notebook

Click through all the cells and run the notebook.

### 3.7 Remember to hit the blue Save button

If you hit the blue Save button it saves the machine image for you and you can stop your notebook and resume it when you are ready and your environment will be saved.

Remember to use github to save your changes.  Work on a branch and ask for pull requests.
Or be sure you are not working on a figure that someone else is working on so we don't have conflicts.  Ask if you have questions!

![click sbas folder](https://github.com/TheJacksonLaboratory/sbas/blob/master/assets/RememberSaveJupyterSessionWithBlueSaveButtonAndSaveOnGithubYourChanges.png)

As you click, the data files are unpacked
### 
## 4 Content Details of files

You can get the data from the jupyternotebooks if you use piggyback or we have the in the following buckets, in the cloud environment.
To access the results from differential gene expression analysis and differential splicing junction analysis.


The results from rmats-nf NextFlow pipeline analysis of the GTEx V8 release are found in github release repositories but for convenience of our analysis now are in the google bucket named robinson-bucket in the folder named notebooks, organised in .tar.gz files


There are 5 files.

| filename | Description  |
| ------------ | ------------- |
| data.tar.gz | Contains results from the rmats-nf pipeline and the results from the execution of the first differential gene expression pipeline as captured in sbas/jupyter/differentialGeneExpressionAnalysis.ipynb |
| assets.tar.gz | This file contains the tissues.tsv file that was used for generating display labels (see the code for how to do this in the current sbas/jupyter/figure1.ipynb notebook - soon to be named sexBiasedASEventsPerTissue.ipynb). This reduces our tissue list to 39 based upon number of samples with adequate data.  Selection is indicated by the column `include` where 1 indicates keep and 0 indicates exclude |
| DGE_splicing_data.tar.gz | Results from the execution of the differential splicing junction analysis as executed for each of the reduced 39 tissue as outlined in assets.tar.gz.   The file details are explained in the next section |
| DGE_splicing_pdf.tar.gz | These are the MDS_plots -- useful for investigating the ability of voom to clean up the segregation of results -- a summary figure will be included in the publication of some subset of this information. |

## 3.1 Contents of data.tar.gz

Also available through github using wget or piggyback within R - location noted in box


| filename | Description  |
| ------------ | ------------- |
| DGE_gene_csv.2020-06-17.tar.gz | Results from jupyter/differentialGeneExpressionAnalysis.ipynb |
| fromGTF.tar.gz | Annotation files from rMATS 3.2.5 on the gencode v.30 GTF file |
| gtex.corrected.rds | Corrected Yarn generated expressionSet object on the gene expression data from GTEx v8 with the phenotype data available - |
| rmats_final.a3ss.jc.ijc.txt.gz | Inclusion junction counts for the a3ss |
| rmats_final.a3ss.jc.inclen.txt.gz | rMATS percent spliced in for a3ss |
| rmats_final.a3ss.jc.sjc.txt.gz | Skipped junction counts for a3ss |
| rmats_final.a3ss.jc.skiplen.txt.gz | Length of junction for the alternative a3ss |
| rmats_final.a5ss.jc.ijc.txt.gz | Inclusion junction counts for the a5ss  |
| rmats_final.a5ss.jc.inclen.txt.gz | Inclusion junction length for a5ss |
| rmats_final.a5ss.jc.inc.txt.gz | rMATS percent spliced in for a5ss |
| rmats_final.a5ss.jc.sjc.txt.gz | Skipped junction counts for a5ss |
| rmats_final.a5ss.jc.skiplen.txt.gz | Skipped junction length for a5ss |
| rmats_final.mxe.jc.ijc.txt.gz | Inclusion junction counts for mxe |
| rmats_final.mxe.jc.inclen.txt.gz | length of junction for the mxe |
| rmats_final.mxe.jc.inc.txt.gz | rMATS percent spliced in for mxe |
| rmats_final.mxe.jc.sjc.txt.gz | Skipped junction counts mxe |
| rmats_final.mxe.jc.skiplen.txt.gz | Skipped junction length for mxe |
| rmats_final.ri.jc.ijc.txt.gz | Inclusion junction counts for ri |
| rmats_final.ri.jc.inclen.txt.gz | Inclusion junction length for ri |
| rmats_final.ri.jc.inc.txt.gz | rMATS percent spliced in for ri |
| rmats_final.ri.jc.sjc.txt.gz |  Skipped junction counts for ri |
| rmats_final.ri.jc.skiplen.txt.gz | Skipped junction length |
| rmats_final.se.jc.ijc.txt.gz | Inclusion junction counts for se |
| rmats_final.se.jc.inclen.txt.gz | Inclusion junction length for se |
| rmats_final.se.jc.inc.txt.gz | rMATS percent spliced in for se |
| rmats_final.se.jc.sjc.txt.gz | Skipped junction counts for se |
| rmats_final.se.jc.skiplen.txt.gz | Skipped junction length |
| SraRunTable.txt.gz | Annotation file connecting SRR with SAMPID |
| srr_pdata.csv | Updated phenodata matching SAMPID to SRR (1 to many) |
| Tissue_list.csv | Listing of tissues for this study in assets put as part of this release for convenience |
