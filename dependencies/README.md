# Create a requirements.txt from the output of `conda list`

1. Capture conda installed packages 

```bash 
conda list > conda_list.txt
```

2. Load in R the dataframe with installed packages metadata

```r
library(readr)
conda_list <- read_table2("dependencies/conda_list.txt"
```

3. Convert metadata to required syntax formatted

Following the convention for conda packages:

`channel::package=version`

```r
conda_list$requirements<- paste0(conda_list$Channel, "::", conda_list$Name, "=", conda_list$Version)
```

4. Write the `requirements` column as a txt file:

```r
library(data.table)
data.table::fwrite(conda_list[,4], file = "dependencies/requirements.txt", quote = FALSE)
```
