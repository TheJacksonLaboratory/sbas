Sys.setenv(TAR = "/bin/tar")
BiocManager::install(c('bumphunter' , 'GenomeInfoDb', 'GenomeInfoDbData')
install.packages(c('doParallel', 'doRNG', 'foreach', 'stringi'), repo = 'https://cran.r-project.org')
devtools::install_github("TheJacksonLaboratory/yarn", ref = "annes-changes", upgrade="never")