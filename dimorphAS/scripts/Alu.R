
setwd('/Users/karleg/Dimorph/')

an.tab<-read.table('/Users/karleg/Dimorph/fromGTF.SE.txt',sep='\t',header=T)

up.intron.s<-an.tab$upstreamEE

up.intron.e<-an.tab$downstreamES

out.tab<-cbind(as.character(an.tab$chr),up.intron.s,up.intron.e,as.character(an.tab$geneSymbol),up.intron.e-up.intron.s-1,as.character(an.tab$strand),an.tab$ID)

write.table(out.tab,'all_upstream_introns.bed',sep='\t',col.names = F,row.names = F,quote = F)

system('bedtools intersect -wa -wb -s -F 1 -a all_upstream_introns.bed -b repeats.bed > all_repeats_intersections.bed')

alu.tab<-read.table('/Users/karleg/Dimorph/all_repeats_intersections.bed',sep='\t',header = F)

colnames(alu.tab)<-c('chr','up_intr_s','up_intr_e','gene','up_intr_l','strand','ID','chr','alu_s','alu_e','alu','blank','strand')

alu.tab<-alu.tab[grepl('^Alu',alu.tab$alu),]

alu.tab<-alu.tab[!duplicated(paste0(alu.tab$ID,alu.tab$alu)),]

all.tab<-table(alu.tab$alu)

length(all.tab)

total.as<-nrow(an.tab)

all.sb.events<-c()

all.events<-c()

meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.sets<-list(c("Brain - Cerebellum","Brain - Cerebellar Hemisphere"),
                  c("Brain - Caudate (basal ganglia)","Brain - Nucleus accumbens (basal ganglia)","Brain - Putamen (basal ganglia)"),
                  c("Brain - Cortex","Brain - Frontal Cortex (BA9)","Brain - Anterior cingulate cortex (BA24)"),
                  c("Brain - Amygdala","Brain - Hippocampus"),
                  c("Esophagus - Gastroesophageal Junction","Esophagus - Muscularis"),
                  c("Skin - Sun Exposed (Lower leg)","Skin - Not Sun Exposed (Suprapubic)")
)
for (singleton.set in unique(as.character(meta.data$body_site_s[!meta.data$body_site_s %in% unlist(tissue.sets)])))
  
  tissue.sets[length(tissue.sets)+1]<-singleton.set

for (tissue.set in tissue.sets)
{
  if (file.exists(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),'se.txt')))
    
    sb.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),'se.txt'))
  
  all.sb.events<-unique(c(all.sb.events,rownames(sb.tab)[abs(sb.tab$logFC)>=(log2(1.5)) & sb.tab$adj.P.Val<=0.05]))
  
  all.events<-unique(c(all.events,rownames(sb.tab)))
}

all.tab<-table(alu.tab$alu[alu.tab$ID %in% all.events])

length(all.tab)

total.as<-length(all.events)

total.sb<-length(all.sb.events)

sb.tab<-table(alu.tab$alu[alu.tab$ID %in% all.sb.events])

res<-matrix(ncol=7,nrow=0)

colnames(res)<-c('Alu','CountAS','CountSB','TotalAS','TotalSB','PVal','FDR')

for (i in (1:length(all.tab)))
{
  
  cur.cat<-names(all.tab)[i]  
  
  if (!cur.cat %in% names(sb.tab) | all.tab[cur.cat]<10 | !grepl('^Alu',cur.cat))
    
    next
  
  count.sb<-sb.tab[cur.cat]
  
  count.as<-all.tab[cur.cat]
  
  pval<-exp(phyper(q=count.sb, m=count.as, n=total.as-count.as, k=total.sb,log.p = T,lower.tail = F))
  
  res<-rbind(res, c(cur.cat,count.as,count.sb,total.as,total.sb,pval,0))
  
  
}

res[,ncol(res)]<-mt.rawp2adjp(as.numeric(res[,'PVal']),proc='BH')[[1]][,2][order(mt.rawp2adjp(as.numeric(res[,'PVal']),proc='BH')$index)]

View(res[order(as.numeric(res[,'PVal'])),])


alu.sum<-read.table('/Users/karleg/Downloads/alu_summary_robin.txt',sep='\t',header = F)

alu.sum[,2]<-gsub('.bed','',alu.sum[,2])

merged.tab<-merge(alu.sum,res,by.x='V2',by.y='Alu')

for (i in c(8,11,14,17,19))
  
  merged.tab[,i]<-as.numeric(merged.tab[,i])

plot(-log10(as.numeric(as.character(merged.tab$PVal))), rowMeans(merged.tab[,c(8,11,14,17,19)]),xlab='-log10(P-val)',ylab='Mean Z')

y=rowMeans(merged.tab[,c(8,11,14,17,19)]) 
  
x=(-log10(as.numeric(as.character(merged.tab$PVal))))

fit<-lm(y ~ x)

abline(fit)

summary(fit)

write.table(res[order(as.numeric(res[,'PVal'])),],'alu_enrichment.txt',sep='\t',row.names = F,col.names = T,quote = F)




esr.tab<-read.table('/Users/karleg/Downloads/all-byclass-noM-plusminus-hits-uniq.txt',sep=' ',header = T)

chr<-unlist(lapply(lapply(lapply(lapply(as.character(esr.tab$item_coordinates),strsplit,split=':'),'[[',1),'[[',1),'[',1))

start<-unlist(lapply(lapply(lapply(as.character(esr.tab$item_coordinates),strsplit,split=':'),'[[',1),'[',2))

start<-unlist(lapply(lapply(lapply(start,strsplit,split='-'),'[[',1),'[',1))

end<-unlist(lapply(lapply(lapply(as.character(esr.tab$item_coordinates),strsplit,split=':'),'[[',1),'[',2))

end<-unlist(lapply(lapply(lapply(end,strsplit,split='-'),'[[',1),'[',2))

end<-unlist(lapply(lapply(lapply(end,strsplit,split=','),'[[',1),'[',1))

strand<-unlist(lapply(lapply(lapply(as.character(esr.tab$item_coordinates),strsplit,split=':'),'[[',1),'[',2))

strand<-unlist(lapply(lapply(lapply(strand,strsplit,split=','),'[[',1),'[',2))

out.tab<-cbind(chr,start,end,rep('NA',length(chr)),as.integer(end)-as.integer(start)-1,strand,as.character(esr.tab$sequence_name),esr.tab$item_antisense)

write.table(out.tab,'alu_with_esr.bed',sep='\t',col.names = F,row.names = F,quote = F)

system('bedtools intersect -wa -wb -s -F 1 -a all_upstream_introns.bed -b alu_with_esr.bed > all_esr_intersections.bed')

rm(out.tab)

rm(esr.tab)

rm(chr)

rm(start)

rm(end)

rm(strand)

alu.esr.tab<-read.table('all_esr_intersections.bed',sep='\t',header = F)

alu.esr.tab<-alu.esr.tab[grepl('^Alu',alu.esr.tab$V14),]

alu.esr.tab<-alu.esr.tab[TRUE==alu.esr.tab$V15,]

alu.esr.tab<-alu.esr.tab[!duplicated(paste0(alu.esr.tab$V7,alu.esr.tab$V14)),]

alu.esr.tab<-alu.esr.tab[alu.esr.tab$V7 %in% all.events,]

alu.types<-table(as.character(alu.esr.tab$V14))

sb.esr.tab<-table(as.character(alu.esr.tab$V14[alu.esr.tab$V7 %in% all.sb.events]))

res2<-matrix(ncol=7,nrow=0)

colnames(res2)<-c('Alu','CountAS','CountSB','TotalAS','TotalSB','PVal','FDR')

for (i in (1:length(alu.types)))
{
  
  cur.cat<-names(alu.types)[i]  
  
  count.sb<-sb.esr.tab[cur.cat]
  
  count.as<-alu.types[cur.cat]
  
  if (is.na(count.sb) || all.tab[cur.cat]+sb.tab[cur.cat]<25 || count.sb<25)
    
    next
  
  pval<-exp(phyper(q=count.sb, m=count.as, n=all.tab[cur.cat]-count.as, k=sb.tab[cur.cat],log.p = T,lower.tail = F))
  
  res2<-rbind(res2, c(cur.cat,count.as,count.sb,all.tab[cur.cat],sb.tab[cur.cat],pval,0))
}

res2[,ncol(res2)]<-mt.rawp2adjp(as.numeric(res2[,'PVal']),proc='BH')[[1]][,2][order(mt.rawp2adjp(as.numeric(res2[,'PVal']),proc='BH')$index)]

View(res2[order(as.numeric(res2[,'PVal'])),])

write.table(res2[order(as.numeric(res2[,'PVal'])),],'esr_enrichment.txt',sep='\t',quote = F,row.names = F)

