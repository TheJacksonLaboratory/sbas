#!/usr/bin/env Rscript

############################## ARGUMENTS SECTION #############################
## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no all arguments passed or help needed
if("--help" %in% args | (length(args) == 0) ) {
  cat("
      The R Script differential_expression.R

      Mandatory arguments:
          --counts=type             - description
          --metadata=type           - description
          --help                    - you are reading it

      Optionnal arguments:
          --threshold=Value         - a threshold, default:0.25

      Usage:
      
          The typical command for running the script is as follows:
    
          ./differential_expression.R --counts='testdata/counts.csv' --metadata='testdata/metadata.csv' --threshold=0.20 
      
      
      WARNING : here put all the things the user has to know
      \n")

  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))
names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)

## Give some value to optional arguments if not provided
if(is.null(args$threshold)) {args$threshold =0.25} else {args$threshold=as.numeric(args$threshold)}

cat("\n")
cat("counts     : ", args$counts,"\n",sep="")
cat("metadata   : ", args$metadata,"\n",sep="")
cat("threshold  : ", args$threshold,"\n",sep="")


counts_df <- readr::read_csv(file = args$counts)
counts_df[1:2, 1:2]
