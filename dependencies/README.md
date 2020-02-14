# Create a `requirements.txt` from the output of `conda list -e`

1. Capture conda installed packages 

```bash
conda list > dependencies/conda_list.txt
```    

2. Reformat to create a requirement.txt files

```R 
library(readr)
conda_list              <- read_table2("dependencies/conda_list.txt", col_names = FALSE, comment = "#")
colnames(conda_list)    <- c('package', 'version', 'build', 'channel')
conda_list$requirements <- paste0(conda_list$channel, "::", conda_list$package, "=", conda_list$version)
data.table::fwrite(conda_list[,'requirements'], file = "dependencies/requirements.txt", quote = FALSE, col.names = FALSE)
```


# Update the `base` conda environment

```bash

conda env update --name base --file dependencies/requirements.txt 

```

#  Install `TheJacksonLaboratory/yarn`


```bash

Rscript dependencies/install.R 

```
