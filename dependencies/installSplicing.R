Sys.setenv(TAR = "/bin/tar")
install.packages(c('BiocManager', 'snakecase','multtest','doParallel', 'readr', 'doRNG', 'foreach', 'stringi', 'pheatmap', 'R.utils'), repo = 'https://cran.r-project.org')
BiocManager::install(c( 'minfi','bumphunter' , 'GenomeInfoDb', 'GenomeInfoDbData', 'edgeR', 'statmod'))
