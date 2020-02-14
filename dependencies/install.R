Sys.setenv(TAR = "/bin/tar")
install.packages(c('doParallel', 'doRNG', 'foreach'), repo = 'https://cran.r-project.org')
devtools::install_github("TheJacksonLaboratory/yarn", ref = "annes-changes", upgrade="never")
