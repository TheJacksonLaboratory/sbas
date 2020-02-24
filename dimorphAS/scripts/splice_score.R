library(stringr)

library(ggplot2)

my.grepl<-function(x)
{
  
  grepl(as.character(unlist(x['ACCESSION_ID'])),as.character(unlist(x['nmd.ids'])))
  
}

s.tab<-read.csv('/Users/karleg/Downloads/SplicingCalculatorData/hg38_ensembl_splicing_calculator_output.tsv',sep='\t',header = T)

s.tab<-s.tab[s.tab$DONOR!='NaN' & s.tab$ACCEPTOR!='NaN',]

##keep only exons that are in both tables

s.tab<-cbind(s.tab,paste(paste0('chr',s.tab$X.CHR),s.tab$BEGIN,s.tab$END,s.tab$STRAND,sep='-'))

colnames(s.tab)[ncol(s.tab)]<-'EXON_ID'

nmd.tab<-read.table('/Users/karleg/Dimorph/NMD_summary.txt',sep='\t')

an.tab<-read.table('/Users/karleg/Dimorph/fromGTF.SE.txt',header=T)

nmd.tab<-cbind(nmd.tab,paste(an.tab$chr,an.tab$exonStart_0base,an.tab$exonEnd,an.tab$strand,sep='-'))

colnames(nmd.tab)[ncol(nmd.tab)]<-'EXON_ID'

nmd.tab<-nmd.tab[nmd.tab$num.transcripts>0,]

s.tab<-s.tab[s.tab$EXON_ID %in% nmd.tab$EXON_ID,]

s.tab<-merge(s.tab,nmd.tab,by.x='EXON_ID',by.y='EXON_ID')

s.tab<-s.tab[unlist(lapply(split(s.tab[,c('ACCESSION_ID','nmd.ids')],1:nrow(s.tab)),my.grepl)) | s.tab$num.nmd==0,]

s.tab<-s.tab[!duplicated(paste(s.tab$EXON_ID,s.tab$num.nmd)),]

mean((s.tab$ACCEPTOR+s.tab$DONOR)[s.tab$num.nmd>0])

mean((s.tab$ACCEPTOR+s.tab$DONOR)[s.tab$num.nmd==0])

df<-data.frame(type=ifelse(s.tab$num.nmd>0,'NMD','Non-NMD'),score=s.tab$ACCEPTOR+s.tab$DONOR)

p<-ggplot(df, aes(x=type, y=score,color=type))+ geom_boxplot()+ylab('Acceptor+Donor Scores')+xlab('')+theme(axis.text=element_text(size=16))+theme(plot.background = element_rect(size=0.5,linetype="solid",color="black"))

plot(p)

wilcox.test((s.tab$ACCEPTOR+s.tab$DONOR)[s.tab$num.nmd==0],(s.tab$ACCEPTOR+s.tab$DONOR)[s.tab$num.nmd>0])

