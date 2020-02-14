# Create a `requirements.txt` from the output of `conda list -e`

1. Capture conda installed packages 

This file may be used to create an environment using:
conda create --name <env> --file <this file>

```bash 
conda list -e > dependencies/requirements.txt 
```
