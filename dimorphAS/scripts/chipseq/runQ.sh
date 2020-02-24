#!/bin/bash

###  Script to download a ChIP-seq data set from ###
###  NIH SRA, align it, and do Q ChIP-seq peak   ###
###  calling.                                    ###


###  Section 1. Required Programs/Paths          ###
## obviously, adjust paths as needed

# 1. fastq-dump
# download here: https://ncbi.github.io/sra-tools/install_config.html
#tar xvfz ncbi-sra-tools-2.8.2-4-0-g0597af0.tar.gz (source) or
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software (precompiled for linux--recommended)
# tar xvfz /home/robinp/bin/sratoolkit.2.8.2-1-ubuntu64.tar.gz
# The fastq-dump binary is in the bin subdirectory after unpacking
FASTQDUMP=${HOME}/bin/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump

## We will use bowtie2 for the alignment. bwa would work just as well
## Bowtie2 can be installed on debina/ubuntu using apt-get install bowtie2. This is what I did
BOWTIE2=/usr/bin/bowtie2
## We also need a genome index to run bowtie2.
# We can make the index ourselves or we can use the index from the bowtie2 website
# http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
# For this example, I will use the human hg19 index because I already have it locally
# for real, we should use the hg38 genome, but everything should be completely analogous
#Note when you unpack the index file it has lots of files, hg19.1.bt2  hg19.2.bt2  hg19.3.bt2  hg19.4.bt2  hg19.rev.1.bt2  hg19.rev.2.bt2
# we refer to simply hg19 to get all of these files, as follows (I have all the files in bin/bowtie2):
BT2INDEX=${HOME}/bin/bowtie2/hg19
## I will also assume that samtools is available; it can also be installed using apt-get
SAMTOOLS=/usr/bin/samtools
## To download Q, use the GitHub repository
##  git clone https://github.com/charite/Q.git
# then enter $./configure; make
## Here is the path on my system
Qexecutable=${HOME}/GIT/Q/bin/Q


###  Section 2. Download ChIP-seq dataset from SRA ###
## NOTE-- either use fastqdump to download or go to the SRA FTP
## fastq-dump can be made automatic but it is not reliable
## Explore here
## https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/
## This can certainly be scripted given the SRR accession numbers
## although I have not tried to script getting hundreds of files from here yet

# We want to get to files from the SRA study
# ChIP-seq Analysis of MCF7 Cells on Treatment with Estradiol with internal controls:
# https://www.ncbi.nlm.nih.gov/sra/SRX3447358[accn]
# 1. SRR6350627	// 1 run, 30.3M spots, 1.5G bases, 531.3Mb  SLX-15091_Input_45Min
# 2. SRR6350628	// 31,140,838	1.6G	547.9Mb	2017-12-15  SLX-15091_1b_ER_45min; Homo sapiens; ChIP-Seq

# So we will use these two as a paired sample, one with the specific anti-ER antibody
## and one input (no specific AB, i.e., control)

# Since both files start with SRR635 I went to that part of the SRA acrhive
# https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR635/
# There is one directory for each file, e.g.,
# https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR635/SRR6350627/
# contains SRR6350627.sra
# Therefore, the following two calls to wget


if [ -e SRR6350627.sra ]; then
    echo "Cowardly refusing to re-download SRR6350627.sra"
else
    wget --progress=bar https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR635/SRR6350627/SRR6350627.sra
fi

if [ -e SRR6350628.sra ]; then
    echo "Cowardly refusing to re-download SRR6350628.sra"
else
    wget --progress=bar https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR635/SRR6350628/SRR6350628.sra
fi

for sample in SRR6350627.sra SRR6350628.sra
do
    FASTQ="${sample/.sra/_pass_1.fastq.gz}"
    if [ -e $FASTQ ]; then
	echo "Cowardly refusing to repeat fastqdump extraction"
    else
	echo "Running fastq-dump on sample $sample (may take tens of minutes)"
	$FASTQDUMP --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip $sample
    fi
done

###  Section 3. Align ChIP-seq data with bowtie2  ###
## From section 2, we have two single-end samples in fastq.gz format. We need to align them to the human genome
## Note that bowtie2 can run on multiple threads (--threads <int>),
# see the WebPage for details! I let it run on 1 thread here.

for FASTQ in SRR6350627_pass_1.fastq.gz SRR6350628_pass_1.fastq.gz
do
    OUTNAME="${FASTQ/_pass_1.fastq.gz/.sam}"
    if [ -e ${OUTNAME} ]; then
	echo "Cowardly refusing to re-create alignment file ${OUTNAME}"
    else
	echo "Will run bowtie2 to create alignment file $OUTNAME"
	${BOWTIE2} -x ${BT2INDEX} -U ${FASTQ} -S ${OUTNAME}
    fi
done


### Section 4. Convert SAM to BAM and sort ####
### We run the ChIP-seq peak called with a sorted BAM file, which we create in this step


for SAM in SRR6350627.sam SRR6350628.sam
do
    BAM="${SAM/.sam/_sorted}" ## Note samtools adds suffix automatically
    BAM_WITH_SUFFIX="${SAM/.sam/_sorted.bam}"
    if [ -e ${BAM_WITH_SUFFIX} ]; then
	echo "cowardly refusing to re-create BAM file $BAM"
    else
	echo "Creating BAM file $BAM..."
	${SAMTOOLS} view -bS ${SAM} | samtools sort - ${BAM}
    fi
done


### Section 5. Run Q on the two BAM files  ####
## Note that there are experimental (with specific antibodies) and input (control) files for ChIP-seq
## Any larger scale pipeline will need to keep track of pairs of files for experiment and control files.
## I am doing this manually here.
## There is no way I know of that we can do this automatically with SRA, so we will have to agree on
## a curation format


echo "We will now run Q to call ChIP-seq peaks on the sorted BAM data ..."

${Qexecutable} -t SRR6350628_sorted.bam -c SRR6350627_sorted.bam -o SRR635062_7_8called -wbt -wbc

    



