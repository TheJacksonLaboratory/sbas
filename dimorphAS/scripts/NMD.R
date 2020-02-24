library(Biostrings)

event.type='SE'

setwd('/Users/karleg/Dimorph/')

fasta.file='/Users/karleg/STAR/STAR/bin/MacOSX_x86_64/data/GRCh38/sequence/GRCh38_r91.all.fa'

from.gtf<-read.table(paste0('/Users/karleg/Dimorph/fromGTF.',event.type,'.txt'),header=T)

res<-matrix(c(rep(0,nrow(from.gtf)),rep(0,nrow(from.gtf)),rep('',nrow(from.gtf))),ncol=3)

colnames(res)<-c('num.nmd','num.transcripts','nmd.ids')

finished=0

for (chr in as.character(unique(from.gtf$chr)))
{
  
  if (chr=='chrY')
    
    next
  
  system(paste0('grep ^',gsub('chr','',chr),' /Users/karleg/STAR/STAR/bin/MacOSX_x86_64/data/GRCh38/annotation/Homo_sapiens.GRCh38.91.gtf > cur_chrom.gtf'))
  
  cur.gtf<-read.table('cur_chrom.gtf',sep='\t',quote="")
  
  for (exon.itr in ((1:nrow(from.gtf))[from.gtf$chr==chr]))
  {
    exon.rows<-which((cur.gtf$V4==from.gtf$exonStart_0base[exon.itr]+1) & (cur.gtf$V5==from.gtf$exonEnd[exon.itr]) & cur.gtf$V3=='exon')
    
    for (exon.row in exon.rows)
    {
      
      transcript.first.row<-max(which((cur.gtf$V3=='transcript') & ((1:nrow(cur.gtf))<exon.row)  ))
      
      if (sum((cur.gtf$V3 %in% c('transcript','gene')) & ((1:nrow(cur.gtf))>exon.row)  )==0)
      {
        transcript.last.row<-nrow(cur.gtf)
      }else{
      transcript.last.row<-min(which((cur.gtf$V3 %in% c('transcript','gene')) & ((1:nrow(cur.gtf))>exon.row)  ))-1
      }
      
      out.gtf<-cur.gtf[transcript.first.row:transcript.last.row,]
      
      if (sum(c(sum((out.gtf$V4==from.gtf$upstreamES[exon.itr]+1) & out.gtf$V5==from.gtf$upstreamEE[exon.itr])>0, 
      
      sum((out.gtf$V4==from.gtf$downstreamES[exon.itr]+1) & out.gtf$V5==from.gtf$downstreamEE[exon.itr])>0))<2)
        
        next
      
      if (sum(out.gtf$V3=='CDS')<3)
        
        next
      
      write.table(out.gtf,"transcript.gtf",sep='\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
      
      system(paste('./gffread -y gene.fa -g ',fasta.file,' transcript.gtf',sep=''))
      
      seq<-readAAStringSet('gene.fa')
      
      if (length(seq)==0)
        
        next
      
      l.inc<-length(seq[[1]])
      
      out.gtf<-out.gtf[(out.gtf$V4!=(from.gtf$exonStart_0base[exon.itr]+1)) & (out.gtf$V5!=(from.gtf$exonEnd[exon.itr])),]
      
      write.table(out.gtf,"transcript.gtf",sep='\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
      
      system(paste('./gffread -y gene.fa -g ',fasta.file,' transcript.gtf',sep=''))
      
      seq<-readAAStringSet('gene.fa')
      
      l.skip<-length(seq[[1]])
      
      skip.exon.aa.length<-(from.gtf$exonEnd[exon.itr]-(from.gtf$exonStart_0base[exon.itr]+1))/3
        
      res[exon.itr,2]<-as.integer(res[exon.itr,2])+1
      
      if (l.inc<(l.skip+skip.exon.aa.length-1))
      {
        res[exon.itr,1]<-as.integer(res[exon.itr,1])+1
        
        res[exon.itr,3]<-paste(res[exon.itr,3],cur.gtf$V9[transcript.first.row],sep='***')
      
      }
    }
    
  }
  
  system('rm cur_chrom.gtf')
  
  finished<-finished+sum(from.gtf$chr==chr)
  
  print(paste0("Finished: ",finished))
  
}


write.table(res,"NMD_summary.txt",sep='\t',row.names = TRUE,col.names = TRUE,quote = FALSE)