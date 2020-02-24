p = read.csv('/projects/robinson-lab/dimorphAS/results/limma.subset/p_matrix.csv', header = T, row.names = 1)
#rownames(p) = p$Gene
#p <- head(p)

get.roP<-function(p,bp=NULL,rth) {
  k<-ncol(p)
  rnum<-which(apply(p,1,function(x) !any(is.na(x))))
  pval<-stat<-rep(NA,nrow(p))
  if(!is.null(bp)){
    p<-t(apply(p,1, sort,na.last = T))
    bp<-t(apply(bp,1,sort,na.last = T))
    Ug<-matrix(NA,nrow(p),1)
    Ug[rnum,1]<-p[rnum,rth]	
    Ubg<-matrix(bp[,rth],nrow(p),nrow(bp)/nrow(p))
    pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=Ug,pval=pval, FDR=qval)
  }else{
    pval[rnum]<-apply(p[rnum,],1,function(x) {
      ## calculate rOP p-value 
      pbeta(x[order(x)][rth],rth,k-rth+1)
    }) 
    stat[rnum]<-apply(p[rnum,],1,function(x)x[order(x)][rth])			
    qval<-p.adjust(pval,method="BH")
    res<-list(stat=stat,pval=pval,FDR=qval)	
  }
  names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
  return(res)
}

## run rthOP with rth value from 1-45
for (i in 1:45) {
  ofile <- paste(c('/projects/robinson-lab/dimorphAS/results/limma.subset/rthOP/rthOP.', i, '.csv'), sep="", collapse="")
  print (ofile)
  cp <- get.roP(p = p, rth = i)
  cp2 <- merge(p, cp, by=0, all=TRUE)
  write.csv(cp2, file=ofile, row.names = FALSE)
} 
