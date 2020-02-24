import argparse
import os.path
import subprocess
from scipy.stats import binom

# check if file exists. If so, return path
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return arg  # return an open file handle

    
parser=argparse.ArgumentParser()
parser.add_argument("-b","--bam",help="BAM file representing ChIP-seq experiment",
                    type=lambda x: is_valid_file(parser, x),required=True)
parser.add_argument("-d","--bed",help="BED file representing ChIP-seq peaks",
                    type=lambda x: is_valid_file(parser, x),required=True)
parser.add_argument("-f","--fasta",help="FASTA file with ref sequence used for the BAM",
                    type=lambda x: is_valid_file(parser, x),required=True)
parser.add_argument("-s","--samtools",help="Path to samtools executable",default='samtools')
parser.add_argument("-o","--outfile",help="outfile name",default='outfile.txt')

args=parser.parse_args()

bam_path=args.bam
bed_path=args.bed
fasta_path=args.fasta
samtools_path=args.samtools
outfile_name=args.outfile

########################################################################
########################################################################
## Need to convert the bed file into a sorted BED3 for samtools
## The bed file produced by Q is not sorted. The following sort command
## produces the same order as bedtools sort and should be more easily
## portable since bedtools may not be installed.
## Write the sorted file with just the first three columns to a
## file called 'sorted_peaks.bed' This BED3 file is the required format
## for samtools -mpileup (BED6 does not work!)
########################################################################
########################################################################

sort_command=["sort","-k", "1,1", "-k2,2n",bed_path]

file = open('sorted_peaks.bed','w') 

proc=subprocess.Popen(sort_command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,encoding="utf-8" )
while True:
    line = proc.stdout.readline().rstrip()
    if not line:
        break
    ar=line.split()
    print (line)
    if (len(ar)<3):
        continue
    file.write('\t'.join(ar[0:3]))
    file.write('\n')

file.close()

bed_path='sorted_peaks.bed'



## for testing
verbose=False
if verbose:
    print("bam",bam_path)
    print("bed",bed_path)
    print("fasta",fasta_path)
    print("samtools",samtools_path)

########################################################################
########################################################################
## Process one line of the pileup
########################################################################
########################################################################
def process_pileup_line(ar):
    chrom=ar[0]
    pos=ar[1]
    refbase=ar[2]
    depth=int(ar[3])
    pup=ar[4]
    #ignore qualities
    ref_count=0
    alt_a=0
    alt_c=0
    alt_g=0
    alt_t=0
    deletion=0
    ins=0
    N=len(pup)
    i=0
    while i<N:
        c=pup[i]
        i=i+1 #move to next char
        if (c=='.' or c==','):
            ref_count=ref_count+1
        elif (c=='^'):
            i=i+1 #skip the next chararacter, which is the quality
        elif (c=='A' or c=='a'):
            alt_a=alt_a+1
        elif (c=='C' or c=='c'):
            alt_c=alt_c+1
        elif (c=='G' or c=='G'):
            alt_g=alt_g+1
        elif (c=='T' or c=='t'):
            alt_t=alt_t+1
        elif (c=='+'):
            deletion=deletion+1
            #the next character gives the length of the deletion string
            dellen_start=i
            while pup[i].isDigit():
                i=i+1
            delnum = pup[dellen_start,i]
            k=int(delnum)
            i=i+k
        elif (c=='i'):
            ins=ins+1
            #the next character gives the length of the deletion string
            dellen_start=i
            while pup[i].isDigit():
                i=i+1
            delnum = pup[dellen_start,i]
            k=int(delnum)
            i=i+k
    nonref_count=alt_a+alt_c+alt_g+alt_t+deletion+ins
    #print("depth",depth,"nonref",nonref_count)
    return [depth,ref_count,nonref_count,alt_a,alt_c,alt_g,alt_t,deletion,ins]

########################################################################
########################################################################
## Analyze the allelic balance of one line
########################################################################
########################################################################
def analyze_balance(depth,ref_count,nonref_count):
    if ref_count < 2:
        return "B"
    p=float(nonref_count)/float(depth)
    homref=binom.pmf(nonref_count,depth,0.02)
    homalt=binom.pmf(nonref_count,depth,0.98)
    hetbal=binom.pmf(nonref_count,depth,0.5)
    hetunbal=0 #default in case p<60% and p>40%
    if (p<=0.4 and p>=0.05):
        p=max(0.05,p)
        hetunbal=binom.pmf(nonref_count,depth,p)
    elif (p>=0.6 and p<=0.95):
        p=min(0.95,p)
        hetunbal=binom.pmf(nonref_count,depth,p)
    maxhom=max(homref,homalt) #dont care which hom
    maxhet=max(hetbal,hetunbal)
    #print("p",p,"homref",homref,"homalt",homalt,"hetbal",hetbal,"hetunbal",hetunbal)
    if maxhom>maxhet:
        return "HOM"
    elif hetbal>hetunbal:
        return "B"
    else:
        return "U"
        
########################################################################
########################################################################
## Run samtools pileup on the required intervals
########################################################################
########################################################################

command=["samtools", "mpileup","--positions",bed_path,"-f",fasta_path,bam_path]
proc=subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,encoding="utf-8" )

file = open(outfile_name,'w') 

# The first two lines should be skipped
#[mpileup] 1 samples in 1 input files
#<mpileup> Set max per-file depth to 8000
header=proc.stdout.readline() 
header=proc.stdout.readline()

verbose=True

while True:
    line = proc.stdout.readline().rstrip()
    if not line:
        break
    ar=line.split()
    #print (line)
    if (len(ar)<5):
        continue
    try:
        #print(line)
        depth,ref_count,nonref_count,alt_a,alt_c,alt_g,alt_t,deletion,ins=process_pileup_line(ar)
        r=analyze_balance(depth,ref_count,nonref_count)
        #if r=='HOM' and not verbose:
         #   continue
        outstr='\t'.join([ar[0],ar[1],r,str(depth),str(ref_count),str(nonref_count)])
        file.write(outstr)
        file.write('\n')
    except Exception as e:
        print (e)


file.close()
