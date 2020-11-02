Sys.setenv(TAR = "/bin/tar")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.r-project.org')

BiocManager::install(c( 'minfi','bumphunter' , 'GenomeInfoDb', 'GenomeInfoDbData', 'edgeR', 'statmod'))
install.packages(c('doParallel', 'doRNG', 'foreach', 'stringi', 'pheatmap', 'R.utils'), repo = 'https://cran.r-project.org')
install.packages("gprofiler2",  repo = 'https://cran.r-project.org')
devtools::install_github("TheJacksonLaboratory/yarn@3ff72c0")
devtools::install_github("ropensci/piggyback@87f71e8")
gprofiler2