#!/bin/bash

###  Script to download a ChIP-seq data set from ###
###  NIH SRA, align it, and do Q ChIP-seq peak   ###
###  calling.                                    ###


## Here, we will download and process just two samples for now.



###  Section 1. Required Programs/Paths          ###
## obviously, adjust paths as needed

# 1. fastq-dump
# download here: https://ncbi.github.io/sra-tools/install_config.html
#tar xvfz ncbi-sra-tools-2.8.2-4-0-g0597af0.tar.gz (source) or
# https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software (precompiled for linux--recommended)
# tar xvfz /home/robinp/bin/sratoolkit.2.8.2-1-ubuntu64.tar.gz
# The fastq-dump binary is in the bin subdirectory after unpacking
FASTQDUMP=${HOME}/bin/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump

# 2. STAR RNA ALigner https://github.com/alexdobin/STAR
STAR=${HOME}/bin/STAR/bin/Linux_x86_64/STAR

# To run STAR, we need to make an index first.
## The genome index is the location where we have the genome file.
## There are many ways to get this. For this script I am assuming that we have a file called hg38.fa with all of the canonical chromosomes
## Ask Peter if an SOP would be helpful

## NOTE THIS SCRIPT WILL NOT WORK UNLESS YOU GET THE hg38.fa FILE FROM SOMEWHERE ELSE!

GENOME_INDEX_DIRECTORY=hg38_overhang_49
GENOME_FA=${HOME}/data/ucsc/hg38/hg38.fa
HG38_GTF=Homo_sapiens.GRCh38.91.gtf

if [ ! -e ${HG38_GTF} ]; then
    wget ftp://ftp.ensembl.org/pub/release-91/gtf/homo_sapiens/Homo_sapiens.GRCh38.91.gtf.gz
    gunzip ${HG38_GTF}.gz
fi


if [ ! -d ${GENOME_INDEX_DIRECTORY} ]; then
    mkdir ${GENOME_INDEX_DIRECTORY}
else
    echo "Directory ${GENOME_INDEX_DIRECTORY} already exists"
fi


echo "${STAR} --runThreadN 4 --runMode genomeGenerate --genomeDir ${GENOME_INDEX} --genomeFastaFiles ${GENOME_FA} --sjdbGTFfile ${HG38_GTF} --sjdbOverhang 49"
${STAR} --runThreadN 4 --runMode genomeGenerate --genomeDir ${GENOME_INDEX} --genomeFastaFiles ${GENOME_FA} --sjdbGTFfile ${HG38_GTF} --sjdbOverhang 49





###  Section 2. Download RNA-seq data from SRA ###
if [ -e SRR2443180.sra ]; then
    echo "cowardly refusing to download SRR2443180.sra"
else
    wget --progress=bar https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR244/SRR2443180/SRR2443180.sra
fi

if [ -e SRR2443186.sra ]; then
    echo "cowardly refusing to download SRR2443180.sra"
else
    wget --progress=bar https://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR244/SRR2443186/SRR2443186.sra
fi








for sample in SRR2443180.sra SRR2443186.sra
do
    FASTQ="${sample/.sra/_pass_1.fastq.gz}"
    if [ -e $FASTQ ]; then
	echo "Cowardly refusing to repeat fastqdump extraction"
    else
	echo "Running fastq-dump on sample $sample (may take tens of minutes)"
	$FASTQDUMP --split-files --gzip --skip-technical --read-filter pass --origfmt --readids --clip $sample
    fi
done





STAR  --runThreadN 4 --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --genomeDir /projects/karleg/projects/lps/sra/star_files/star_indices_GRCh38.91_overhang_49  --sjdbGTFfile /projects/karleg/projects/lps/sra/star_files/Homo_sapiens.GRCh38.91.gtf  --outFileNamePrefix /projects/robinson-lab/lps/bam/${file}  --readFilesIn /projects/robinson-lab/lps/fastq/${file}_1.fastq 
 
