    1  git clone https://github.com/TheJacksonLaboratory/sbas.git
    2  cd sbas
    3  conda install -c bioconda bioconductor-edger -y
    4  conda install -c conda-forge r-readr -y
    5  conda install -c bioconda bioconductor-multtest -y
    6  conda install -c conda-forge r-r.utils -y 
    7  conda install -c conda-forge r-snakecase -y
    8  conda install -c conda-forge r-statmod -y
    9  conda install -c bioconda bioconductor-stan -y
   10  tar xvzf ../mounted-data/5ee271c1143fa00113c7d138-DGE_splicing_data.tar.gz-5ee271c1143fa00113c7d138/robinson-bucket/notebooks/DGE_splicing_data.tar.gz -C data/
   11  tar xvzf ../mounted-data/5ee271e6143fa00113c7d1b2-data.tar.gz-5ee271e6143fa00113c7d1b2/robinson-bucket/notebooks/data.tar.gz -C data/
   12  more data/tissue_list.csv 
   13  ls -l
   14  cd jupyter/
   15  jupytext --to Rmd differentialSplicingJunctionAnalysis.ipynb 
   16  jupytext --to Rmd countGenesAndEvents.ipynb 
   17  emacs countGenesAndEvents.ipynb 
   18  emacs countGenesAndEvents.Rmd 
   19  jupytext countGenesAndEvents.Rmd countGenesAndEvents.ipynb 
   20  jupytext --to ipynb countGenesAndEvents.Rmd
   21  emacs countGenesAndEvents.Rmd
   22  jupytext --to Rmd countGenesAndEvents.ipynb 
   23  emacs countGenesAndEvents.Rmd
   24  jupytext --to ipynb countGenesAndEvents.Rmd
   25  rm countGenesAndEvents.Rmd
   26  conda install -c bioconda bioconductor-rtracklayer
   27  conda update -n base conda
   28  jupytext --to Rmd
   29  jupytext --to Rmd countGenesAndEvents.ipynb 
   30  mv countGenesAndEvents.Rmd ../Rmd/
   31  git add countGenesAndEvents.ipynb 
   32  git add ../Rmd/countGenesAndEvents.Rmd 
   33  git status
   34  git add ../metadata/countGenesAndEvents_*
   35  git push
   36  git status
   37  git set 
   38  get --help
   39  git --help
   40  git push -u https://github.com/TheJacksonLaboratory/sbas
   41  git status
   42  git add countGenesAndEvents.ipynb 
   43  rm countGenesAndEvents.Rmd~ 
   44  rm differentialSplicingJunctionAnalysis.Rmd 
   45  git status
   46  rm keepme.txt 
   47  git commit -m "corrected and updated countGenesAndEvents to work with downloaded files"
   48  git config --global user.email "adeslat@jax.org"
   49  git config --global user.name "Anne Deslattes Mays, PhD"
   50  git commit -m "corrected and updated countGenesAndEvents to work with downloaded files"
   51  git push
   52  ls -l ../data/*.gz
   53  ls -l
   54  cd ../data
   55  ls -l
   56  ls -l *DGE*gene* | more
   57  ls -l *DGE.csv
   58  date
   59  rm *DGE.csv
   60  ls -l * dge
   61  ls -l dge/
   62  ls -l *DGE*
   63  ls -l *DGE.csv
   64  ls -l *.gz
   65  ls -l dge/
   66  cd ..
   67  cd jupyter/
   68  git add summary_per_tissue_diff_expressed.ipynb 
   69  git commit -m "error in cell 27"
   70  git push
   71  conda install -c r r-reshape2 -y
   72  git pull
   73  conda install -c conda-forge r-reshape
   74  conda install -c conda-forge r-reshape -y
   75  git pull
   76  date
   77  ls -lt  ../data | more
   78  history
   79  history > stepsToLogin.MD
