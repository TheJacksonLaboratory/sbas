shuffle.regs<-function(regs,num.itr=25000)
{
  
  prob<-unlist(lapply(regs,length))/sum(unlist(lapply(regs,length)))
  
  for (itr in (1:num.itr))
  {
    flag=TRUE
    
    while(flag)
    {
        idx1<-sample(length(regs),1,prob=prob)
        
        rbp1<-gsub('\\(.*\\)','',sample(regs[[idx1]],1))
        
        idx2=idx1
        
        while(idx2==idx1)
        {
            idx2<-sample(length(regs),1,prob=prob)
            
            rbp2<-gsub('\\(.*\\)','',sample(regs[[idx2]],1))
        }
        
        if (sum(grepl(paste0('^',rbp1,'\\('),regs[[idx2]]))==0 && sum(grepl(paste0('^',rbp2,'\\('),regs[[idx1]]))==0)
          
        {
          
          regs[[idx1]][which(grepl(paste0('^',rbp1,'\\('),regs[[idx1]]))]<-paste0(rbp2,'()')
          
          regs[[idx2]][which(grepl(paste0('^',rbp2,'\\('),regs[[idx2]]))]<-paste0(rbp1,'()')
          
          flag=FALSE
          
        }
          
          
    }
  
  }
  
  
  return(regs)
  
}


hbm.tab<-read.table('/Users/karleg/Dimorph/summary_hbm.txt',header = T,sep='\t')

hbm.tab<-hbm.tab[hbm.tab$Dimorphic=='Yes',]

tmp<-strsplit(as.character(hbm.tab$Sig..RBPs),split = ',')

#for (i in (1:length(tmp)))
  
 # if (length(tmp[i][[1]])>0)
  
  # tmp[i][[1]]<-tmp[i][[1]][!grepl('0.0',tmp[i][[1]])]

num.regs<-unlist(lapply(tmp,length))





all.genes<-read.table('/Users/karleg/Dimorph/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct',sep='\t',header=T,skip=2,colClasses = c(rep("character", 2), rep("NULL", 11688)))

all.genes<-all.genes[!duplicated(all.genes$Description),]

rbp.names<-unique(gsub('_.*','',list.files('/Users/karleg/Dimorph/RBP_PSSMs/')))

rbp.names<-rbp.names[rbp.names %in% all.genes$Description]


m.comb<-matrix(ncol=length(rbp.names),nrow=length(rbp.names))

rownames(m.comb)<-rbp.names

colnames(m.comb)<-rbp.names

tmp<-as.character(hbm.tab$Sig..RBPs)

tmp<-strsplit(as.character(hbm.tab$Sig..RBPs),split = ',')

tmp<-tmp[lapply(tmp,length)>0]

#for (i in (1:length(tmp)))

# if (length(tmp[i][[1]])>0)

# tmp[i][[1]]<-tmp[i][[1]][!grepl('0.0',tmp[i][[1]])]

num.regs<-unlist(lapply(strsplit(tmp,split = ','),length))

#rbp.exp<-read.table('/Users/karleg/Dimorph/rbp_expression.txt',header=TRUE)

#meta.data<-read.table('/Users/karleg/Dimorph/2017December8GTExRNASeqSRARunTable.txt',sep='\t',header=TRUE)

#colnames(rbp.exp)<-gsub('\\.','-',colnames(rbp.exp))

#tissue.list<-c('Heart - Left Ventricle',
 #              'Breast - Mammary Tissue',
  #             'Brain - Cortex','Brain - Frontal Cortex (BA9)','Brain - Anterior cingulate cortex (BA24)',
   #            'Adrenal Gland'
    #           ,'Adipose - Subcutaneous','Muscle - Skeletal','Thyroid','Cells - Transformed fibroblasts','Artery - Aorta','Skin - Sun Exposed (Lower leg)','Skin - Not Sun Exposed (Suprapubic)')

#meta.data<-meta.data[meta.data$body_site_s %in% tissue.list,]

#rbp.exp<-rbp.exp[,c(1,2,which(colnames(rbp.exp) %in% meta.data$Sample_Name_s))]


for (rbp1.i in 1:(length(rbp.names)-1))
  
  for (rbp2.i in (rbp1.i+1):length(rbp.names))
  {
    
    rbp1<-rbp.names[rbp1.i]
    
    rbp2<-rbp.names[rbp2.i]
    
    obs.pairs<-0   
  
    for (i in (1:length(tmp)))

       if (grepl(rbp1,tmp[i]) && grepl(rbp2,tmp[i]))
       
         obs.pairs<- obs.pairs+1
  
    m.comb[rbp1,rbp2]<-obs.pairs
  
    }



r.m.comb<-matrix(rep(0,length(rbp.names)^2),ncol=length(rbp.names),nrow=length(rbp.names))

rownames(r.m.comb)<-rbp.names

colnames(r.m.comb)<-rbp.names

num.shuffles<-100

for (i in (1:num.shuffles))
{
  
  print(paste0('rand itr: ',i))
  
  tmp<-shuffle.regs(tmp)
  
  
  for (rbp1.i in 1:(length(rbp.names)-1))
    
    for (rbp2.i in (rbp1.i+1):length(rbp.names))
    {
      
      rbp1<-rbp.names[rbp1.i]
      
      rbp2<-rbp.names[rbp2.i]
      
      obs.pairs<-0   
      
      for (j in (1:length(tmp)))
        
        if (grepl(rbp1,tmp[j]) && grepl(rbp2,tmp[j]))
          
          obs.pairs<- obs.pairs+1
        
        r.m.comb[rbp1,rbp2]<-r.m.comb[rbp1,rbp2]+obs.pairs
        
    }
}


r.m.comb[upper.tri(m.comb)]<-r.m.comb[upper.tri(m.comb)]/num.shuffles

m.comb[upper.tri(m.comb)]<-m.comb[upper.tri(m.comb)]/r.m.comb[upper.tri(m.comb)]

hist(m.comb[upper.tri(m.comb)])

library(pheatmap)

diag(m.comb)<-1

m.comb[lower.tri(m.comb)]<-t(m.comb[upper.tri(m.comb)])

pheatmap(m.comb)

write.table(m.comb,'/Users/karleg/Dimorph/RDATA/combinatorial_regulation_matrix.txt',sep='\t',quote = F)

library(TFBSTools)

m.seq<-matrix(rep(0,length(rbp.names)^2),ncol=length(rbp.names),nrow=length(rbp.names))

rownames(m.seq)<-rbp.names

colnames(m.seq)<-rbp.names

pwm.files<-list.files('/Users/karleg/Dimorph/RBP_PSSMs/')

for (rbp1.i in 1:(length(rbp.names)-1))
  
  for (rbp2.i in (rbp1.i+1):length(rbp.names))
  {
    
    rbp1<-rbp.names[rbp1.i]
    
    rbp2<-rbp.names[rbp2.i]
    
    pwm.file1<-pwm.files[which(grepl(paste0(rbp1,'_'),pwm.files))]
    
    pwm.file2<-pwm.files[which(grepl(paste0(rbp2,'_'),pwm.files))]
    
    for (i1 in (1:length(pwm.file1)))
      
      for (i2 in (1:length(pwm.file2)))
      {
    
        pwm1<-t(read.table(paste0('/Users/karleg/Dimorph/RBP_PSSMs/',pwm.file1[i1]),header=T,row.names = 1))
        
        pwm2<-t(read.table(paste0('/Users/karleg/Dimorph/RBP_PSSMs/',pwm.file2[i2]),header=T,row.names = 1))
        
        pwm1<-matrix(as.matrix(pwm1),nrow=4)
        
        rownames(pwm1)<-c("A", "C", "G", "T")
        
        pwm1 <- PWMatrix(profileMatrix=pwm1)
        
        pwm2<-matrix(as.matrix(pwm2),nrow=4)
        
        rownames(pwm2)<-c("A", "C", "G", "T")
        
        pwm2 <- PWMatrix(profileMatrix=pwm2)
        
        m.seq[rbp1,rbp2]<-m.seq[rbp1,rbp2]+PWMSimilarity(pwm1, pwm2, method="Euclidean")
    
    }

    m.seq[rbp1,rbp2]<-m.seq[rbp1,rbp2]/(length(pwm.file1)*length(pwm.file2))
  }


m.seq[lower.tri(m.seq)]<-t(m.seq[upper.tri(m.seq)])

pheatmap(m.seq)

save.image('/Users/karleg/Dimorph/RDATA/combinatorial_regs.RData')
