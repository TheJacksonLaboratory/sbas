# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .sh
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Bash
#     language: bash
#     name: bash
# ---

# # A `>_` bash kernel Jupyter Notebook for testing Ontologizer on the output of HBA-DEALS

# ---
# **NOTE**
#
# In order to work in this document you have to have installed the required dependencies for Ontologizer. To check if you have Ontologizer available in your system please first type the following:
#
#
# ```bash
# java -jar Ontologizer.jar --help
# ```
#
# If Ontologizer is not available in you Jupyter Notebook Session, please follow the instructions and execute the script:
#
# ```bash
# bash install_ontologizer.sh
# ```
#
# ---

java -jar Ontologizer.jar --help

# # Retrieve required data for `Ontologizer`
#
# The `ontologizer_test.tar.gz` has 4 files inside:
#
# | | FILE | DESCRIPTION|
# |--|:---|:---|
# |1|**`universe.txt`** | created by writing all the GeneSymbol entries in the HBA-DEALS results table. For retrieving `GeneSymbol`, the Gene column was splitted by `'_'` into `Geneid` and `GeneSymbol`, eg. `ENSG00000004059.11_ARF5` -> `ENSG00000004059.11`, `ARF5`
# |2|**`gene_set.txt`** | created by writing the `GeneSymbol` entries after applying a filtering criterion (for this test,  I used `P` < 0.05 and `ExpLogFc` > 1.2) |
# |3|**`goa_human.gaf`** | downloaded from here: http://current.geneontology.org/annotations/goa_human.gaf.gz |
# |4|**`go.obo`** | downloaded from here: http://purl.obolibrary.org/obo/go.obo |
#
#
#
# The release with the example data in `ontologizer_test.tar.gz` can be found at:<br>
# https://github.com/cgpu/HBA-DEALS/releases/tag/ontologizer.
#
# To get the url of the `ontologizer_test.tar.gz` file right click and `Copy link adreess` as shown in the gif:
#
# ![](http://g.recordit.co/5IcThtkQ6H.gif)
#
#

# download file and decompress archive in a folder name ontologizer_test
# mv contents of ontologizer_test in current working dir and delete empty folder ontologizer_test
wget https://github.com/cgpu/HBA-DEALS/releases/download/ontologizer/ontologizer_test.tar.gz && \
tar -xvzf ontologizer_test.tar.gz -C . && \
mv ontologizer_test/* . && \
rm -r ontologizer_test

# # Run `Ontologizer` 
# based on [`@karleg's Breast-96-Samples.R (ijc+sjc)`](https://github.com/TheJacksonLaboratory/sbas/commit/c5b1ffcebbbde03057cf85e31e9ae4743103df08#diff-9f1b2a73fd1da7a66fe6b2d603f2f55bR156)

java -jar Ontologizer.jar \
    -g go.obo \
    -a goa_human.gaf \
    -s gene_set.txt \
    -p universe.txt \
    -c Term-For-Term \
    -m Benjamini-Hochberg \
    -n

head -1 anno-gene_set-Term-For-Term-Benjamini-Hochberg.txt 

head -4 table-gene_set-Term-For-Term-Benjamini-Hochberg.txt
