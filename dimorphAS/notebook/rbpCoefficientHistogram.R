## set working directory to source file. 
## Make sure to run the perl script "parseSummary.pl" first
# Removal all variables from workspace
rm(list=ls())
library(ggplot2)
library(ggsci)


dat<-read.table("rbpCoeff.txt",header=FALSE)
colnames(dat) <- c("Coefficient", "status")

dat <- dat[dat$status=='biased',]



p2 = ggplot(dat, aes(x=Coefficient, fill=status)) +
 geom_histogram(colour = "black",  position = "dodge",binwidth=0.1) +
  scale_x_continuous(name = "RBP coefficients",
                     breaks = seq(-3, 3, 1),
                     limits=c(-3.5, 3.5)) +
  theme_bw() + theme(text = element_text(size=20),
           axis.text = element_text(size=20, hjust=0.5),
           axis.title = element_text(size=20),legend.position = "none")  
p2_npg = p2 + scale_fill_npg()
p2_npg 

ggsave(file = "rbpCoefficients.pdf", plot = p2_npg)


## Now for expression

dat<-read.table("exprCoeff.txt",header=FALSE)
colnames(dat) <- c("Coefficient", "status")

dat <- dat[dat$status=='biased',]


p3 = ggplot(dat, aes(x=Coefficient, fill=status)) +
  geom_histogram(colour = "black",  position = "dodge",binwidth=0.1) +
  scale_x_continuous(name = "Expression coefficients",
                     breaks = seq(-3, 3, 1),
                     limits=c(-3.5, 3.5)) +
  theme_bw() + theme(text = element_text(size=20),
                     axis.text = element_text(size=20, hjust=0.5),
                     axis.title = element_text(size=20),
                     ,legend.position = "none")  
p3_npg = p3 + scale_fill_npg()
p3_npg 
ggsave(file = "expressionCoefficients.pdf", plot = p3_npg)

dat<-read.table("sexCoeff.txt",header=FALSE)
colnames(dat) <- c("Coefficient", "status")
dat <- dat[dat$status=='biased',]



p4 = ggplot(dat, aes(x=Coefficient, fill=status)) +
  geom_histogram(colour = "black",  position = "dodge",binwidth=0.1) +
  #geom_histogram(aes(y = ..count..), binwidth = 0.2,   position="identity", alpha=0.9) +
  scale_x_continuous(name = "Sex coefficients",
                     breaks = seq(-3, 3, 1),
                     limits=c(-3.5, 3.5)) +
  theme_bw() + theme(text = element_text(size=20),
                     axis.text = element_text(size=20, hjust=0.5),
                     axis.title = element_text(size=20),
                     legend.title = element_blank(),
                     legend.spacing.x = unit(0.2, 'cm'),
                     legend.position="none")  
p4_npg = p4 + scale_fill_npg()
p4_npg 
ggsave(file = "sexCoefficients.pdf", plot = p4_npg)
