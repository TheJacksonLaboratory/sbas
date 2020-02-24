library(edgeR)

setwd('/Users/karleg/Dimorph/')


system('bedtools intersect -wa -wb -s -v -a repeats.bed -b alu_with_esr.bed > alu_no_esr_intersections.bed')

system('bedtools intersect -wa -wb -s -F 1 -a all_upstream_introns.bed -b alu_no_esr_intersections.bed > alu_without_esr_intersections.bed')

system('rm alu_no_esr_intersections.bed')

system('bedtools intersect -wa -wb -s -v -a all_upstream_introns.bed -b repeats.bed > no_alu.bed')


meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

tissue.sets<-list("Breast - Mammary Tissue")

inc.iso.counts.mem<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.ijc.txt',sep=''),header=TRUE)

skip.iso.counts.mem<-read.csv(paste('/Users/karleg/rMATs-GTEx/rmats_final.se.jc.sjc.txt',sep=''),header=TRUE)

for (tissue.set in tissue.sets)
{
  meta.data<-read.table('2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)
  
  sb.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),'se.txt'))
  
  inc.iso.counts<-inc.iso.counts.mem
  
  inc.iso.counts<-inc.iso.counts[inc.iso.counts$ID %in% rownames(sb.tab)[abs(sb.tab$logFC)>=(log2(1.5)) & sb.tab$adj.P.Val<=0.05],]
  
  meta.data$Sample_Name_s<-gsub('-','\\.',meta.data$Sample_Name_s)
  
  meta.data<-meta.data[meta.data$Sample_Name_s %in% colnames(inc.iso.counts),]
  
  rownames(inc.iso.counts)<-inc.iso.counts$ID
  
  inc.iso.counts<-inc.iso.counts[,-1]
  
  inc.iso.counts<-inc.iso.counts[,order(match(colnames(inc.iso.counts),meta.data$Sample_Name_s))]
  
  inc.iso.counts<-inc.iso.counts[,!grepl('11ILO',meta.data$Sample_Name_s) & (meta.data$body_site_s %in% tissue.set)]
  
  skip.iso.counts<-skip.iso.counts.mem
  
  skip.iso.counts<-skip.iso.counts[skip.iso.counts$ID %in% rownames(sb.tab)[abs(sb.tab$logFC)>=(log2(1.5)) & sb.tab$adj.P.Val<=0.05],]
  
  rownames(skip.iso.counts)<-skip.iso.counts$ID
  
  skip.iso.counts<-skip.iso.counts[,-1]
  
  skip.iso.counts<-skip.iso.counts[,order(match(colnames(skip.iso.counts),meta.data$Sample_Name_s))]
  
  skip.iso.counts<-skip.iso.counts[,!grepl('11ILO',meta.data$Sample_Name_s) & (meta.data$body_site_s %in% tissue.set)]
  
  meta.data<-meta.data[!grepl('11ILO',meta.data$Sample_Name_s),]
  
  meta.data<-meta.data[meta.data$body_site_s %in% tissue.set,]
  
  if (sum(colnames(skip.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data))
    
    quit(1)
  
  if (sum(colnames(inc.iso.counts)==meta.data$Sample_Name_s)!=nrow(meta.data))
    
    quit(1)
  
  exprs.mat <- matrix(data=cbind(as.matrix(skip.iso.counts),as.matrix(inc.iso.counts)),ncol=2*ncol(skip.iso.counts),nrow=nrow(skip.iso.counts),dimnames=list(1:nrow(skip.iso.counts),c(colnames(skip.iso.counts),paste(colnames(skip.iso.counts),'2',sep=''))))
  
  keep.events<-rep(0,nrow(exprs.mat))
  
  groups=paste0(as.character(meta.data$sex_s),c(rep('-skip',ncol(skip.iso.counts)),rep('-inc',ncol(inc.iso.counts))))
  
  groups<-groups[colSums(exprs.mat)>0]
  
  exprs.mat<-exprs.mat[,colSums(exprs.mat)>0]
  
  for (group in names(table(groups)))
  
    keep.events<-keep.events + (rowSums(cpm(exprs.mat[,groups %in% group]) > 1) >= 0.25*min(table(groups)))
  
  print(sum(keep.events>2))
  
  rownames(exprs.mat)<-rownames(inc.iso.counts)
  
  exprs.mat<-exprs.mat[keep.events>2,]
  
  min.depths<-rep(10000000000,nrow(exprs.mat))
  
  for (group in names(table(groups)))
    
    min.depths<-apply(cbind(min.depths,rowSums(cpm(exprs.mat[,groups %in% group])>1)),1,min)
        
  annot.tab<-read.table('fromGTF.SE.txt',header=T)
  
  annot.tab<-annot.tab[annot.tab$ID %in% rownames(exprs.mat),]
  
  annot.tab<-annot.tab[order(match(annot.tab$ID,rownames(exprs.mat))),]
  
  sb.tab<-sb.tab[rownames(sb.tab) %in% rownames(exprs.mat),]
  
  sb.tab<-sb.tab[order(match(rownames(sb.tab),rownames(exprs.mat))),]
             
  out.tab<-cbind(rownames(exprs.mat),as.character(annot.tab$geneSymbol),sb.tab$adj.P.Val,min.depths,annot.tab$chr,annot.tab$exonStart_0base,annot.tab$exonEnd,as.character(annot.tab$strand))
  
  colnames(out.tab)<-c('Event','Gene','FDR','Depth Score','chr','event start','event end','strand')
  
  alu.esr.tab<-read.table('all_esr_intersections.bed',sep='\t',header = F)
  
  alu.esr.tab<-alu.esr.tab[grepl('(^AluSp$|^AluSx$)',alu.esr.tab$V14),]
  
  alu.esr.tab<-alu.esr.tab[TRUE==alu.esr.tab$V15,]
  
  alu.esr.tab<-alu.esr.tab[!duplicated(paste0(alu.esr.tab$V7,alu.esr.tab$V14)),]
  
  alu.esr.tab<-alu.esr.tab[alu.esr.tab$V7 %in% out.tab[,'Event'],]
  
  merged.tab<-merge(alu.esr.tab,out.tab,by.x='V7',by.y='Event')
  
  merged.tab<-merged.tab[order(as.integer(merged.tab[,'Depth Score']),decreasing = T),]
  
  merged.tab<-merged.tab[,c('V7','Gene','FDR','Depth Score','chr','event start','event end','strand','V9','V10','V14')]
  
  colnames(merged.tab)[1]<-'Event'
  
  colnames(merged.tab)[c(9,10,11)]<-c('alu start','alu end','alu type')
  
  write.table(merged.tab,paste0('EventsWithAluESR_antisense_',tissue.set),sep='\t',quote = F,col.names = T,row.names = F)
  
  
  sb.tab<-read.table(paste0('/Users/karleg/Dimorph/other/',paste(tissue.set,collapse='.'),'se.txt'))
  
  sb.tab<-sb.tab[sb.tab$adj.P.Val>0.5,]
  
  annot.tab<-read.table('fromGTF.SE.txt',header=T)
  
  annot.tab<-annot.tab[annot.tab$ID %in% rownames(sb.tab),]
  
  annot.tab<-annot.tab[order(match(annot.tab$ID,rownames(sb.tab))),]
  
  out.tab<-cbind(rownames(sb.tab),as.character(annot.tab$geneSymbol),sb.tab$adj.P.Val,annot.tab$chr,annot.tab$exonStart_0base,annot.tab$exonEnd,as.character(annot.tab$strand))
  
  colnames(out.tab)<-c('Event','Gene','FDR','chr','event start','event end','strand')
  
  alu.esr.tab<-read.table('all_esr_intersections.bed',sep='\t',header = F)
  
  alu.esr.tab<-alu.esr.tab[grepl('(^AluSp$|^AluSx$)',alu.esr.tab$V14),]
  
  alu.esr.tab<-alu.esr.tab[FALSE==alu.esr.tab$V15,]
  
  alu.esr.tab<-alu.esr.tab[!duplicated(paste0(alu.esr.tab$V7,alu.esr.tab$V14)),]
  
  alu.esr.tab<-alu.esr.tab[alu.esr.tab$V7 %in% out.tab[,'Event'],]
  
  merged.tab<-merge(alu.esr.tab,out.tab,by.x='V7',by.y='Event')
  
  merged.tab<-merged.tab[,c('V7','Gene','FDR','chr','event start','event end','strand','V9','V10','V14')]
  
  colnames(merged.tab)[1]<-'Event'
  
  colnames(merged.tab)[c(8,9,10)]<-c('alu start','alu end','alu type')
  
  merged.tab<-merged.tab[order(merged.tab$FDR,decreasing = T),]
  
  write.table(merged.tab[1:100,],paste0('EventsWithAluESR_not_antisense_',tissue.set),sep='\t',quote = F,col.names = T,row.names = F)
  
  
  alu.esr.tab<-read.table('alu_without_esr_intersections.bed',sep='\t',header = F)
  
  alu.esr.tab<-alu.esr.tab[grepl('(^AluSp$|^AluSx$)',alu.esr.tab$V11),]
  
  alu.esr.tab<-alu.esr.tab[!duplicated(paste0(alu.esr.tab$V7,alu.esr.tab$V11)),]
  
  alu.esr.tab<-alu.esr.tab[alu.esr.tab$V7 %in% out.tab[,'Event'],]
  
  merged.tab<-merge(alu.esr.tab,out.tab,by.x='V7',by.y='Event')
  
  merged.tab<-merged.tab[,c('V7','Gene','FDR','chr','event start','event end','strand','V9','V10','V11')]
  
  colnames(merged.tab)[1]<-'Event'
  
  colnames(merged.tab)[c(8,9,10)]<-c('alu start','alu end','alu type')
  
  merged.tab<-merged.tab[order(merged.tab$FDR,decreasing = T),]
  
  write.table(merged.tab[1:100,],paste0('EventsWithAluNoESR',tissue.set),sep='\t',quote = F,col.names = T,row.names = F)
  
  
  no.alu.tab<-read.table('no_alu.bed',sep='\t',header = F)
  
  no.alu.tab<-no.alu.tab[no.alu.tab$V7 %in% out.tab[,'Event'],]
  
  merged.tab<-merge(no.alu.tab,out.tab,by.x='V7',by.y='Event')
  
  merged.tab<-merged.tab[,c('V7','Gene','FDR','chr','event start','event end','strand')]
  
  colnames(merged.tab)[1]<-'Event'
  
  merged.tab<-merged.tab[order(merged.tab$FDR,decreasing = T),]
  
  write.table(merged.tab[1:100,],paste0('EventsNoAlu',tissue.set),sep='\t',quote = F,col.names = T,row.names = F)
  
  
  
  rm(exprs.mat)
  
  
}