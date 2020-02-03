#!/bin/sh
inputfile=$1
lines=$(wc -l $1)

print $lines

awk -F ',' '{ for (i=1;i<=NF;i++) sum[i]+=$i} END{for (i in sum) printf("%d ", sum[i]); printf("\n") } $1
