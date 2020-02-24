#This script is use in FindTargetsFromPWMDir.R in order to find RBP binding targets

library(TFBSTools)

library(DBI)

library(Biostrings)

#path to pfm

run.tfbstools<-function(params)
{
#read in the PWM matrix:
pwm.tab<-t(read.table(params[1],header=T,row.names = 1))

#extract the coordinates of the target exon:

chr<-strsplit(params[2],':')[[1]][1]

start<-as.numeric(strsplit(strsplit(params[2],':')[[1]][2],'-')[[1]][1])

end<-as.numeric(strsplit(strsplit(params[2],':')[[1]][2],'-')[[1]][2])

pwm.matrix<-matrix(as.matrix(pwm.tab),nrow=4)

rownames(pwm.matrix)<-c("A", "C", "G", "T")

#create a PWMatrix object: 
pwm <- PWMatrix(profileMatrix=pwm.matrix)

#extract the sequence of the exon +-100bp:
system(paste('samtools faidx /projects/karleg/projects/lps/sra/star_files/GRCh38_r91.all.fa ', paste(chr,':',as.numeric(start)-100,'-',as.numeric(end)-100,sep=''),paste('>exon_seq',chr,':',start,end,'.fa',sep=''),sep=''),intern=TRUE)

#read the exon sequence:
gene.seq<-readDNAStringSet(paste('exon_seq',chr,':',start,end,'.fa',sep=''))
  
subject <- DNAString(gene.seq[[1]])

#predict binding sites
siteset <- searchSeq(pwm, subject, min.score="80%", strand="*")

file.remove(paste('exon_seq',chr,':',start,end,'.fa',sep=''))

return(siteset)

}
#write.table(cbind(writeGFF3(siteset, scoreType="absolute"),pvalues(siteset)),paste('targets/event',params[3],basename(params[1]),'sites.txt',sep='_'),sep='\t',quote = FALSE,row.names = FALSE,col.names = TRUE)




