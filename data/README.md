## The folder `data`

We use the folder **`data`** to store the artefacts that are generated during the analysis, so this folder should appear empty on GitHub. Large files cannot be pushed to GitHub as done with code. One different approach we can adopt for tracking medium sized artefacts up to 2GB. For this purpose we utilise the R package [`ropensci/piggyback`](https://github.com/ropensci/piggyback). This is still a work in progress and we will only sustain this for developing for as long as parts of this projects are a work in progress. See the following section about utilising cached intermediate files for more details.

sha256sum gtex.corrected.rds

c1d5da42c902d7d43e275f59c0d119933404ce43d32e6388696061c5faeede1d  gtex.corrected.rds