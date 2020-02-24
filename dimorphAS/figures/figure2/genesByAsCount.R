## Create Figure 2d -- genes with >10 AS events in multiple tissues
# Removal all variables from workspace
rm(list=ls())

library(ggplot2)
library(dpylr)
library(magrittr)


#ADDED
dat<-read.table("genesWithCommonAs.tsv",header=TRUE)
res <- dat[dat$n>10,]  #remove rows with less than 10 events
res %>% arrange(desc(n)) # rearrange descending according to number of affected tissues
res$GeneSymbol <- factor(res$GeneSymbol, levels = res$GeneSymbol)



# Genes with more than 10 splicing events
p<-ggplot(res, aes(x = GeneSymbol, y = n)) +
  geom_point(size = 4, aes(fill = Tissues, color = Tissues)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=18, angle = 270, hjust = 0.0, vjust = 0.5),
        axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face="plain", colour="black", 
                                    size=24),
        legend.title=element_blank(),
        legend.text = element_text(face="plain", colour="black", 
                                   size=18)) +
  scale_fill_viridis_c(aesthetics = c("colour", "fill"),
                       option = "plasma",
                       limits = c(1, 30), breaks = c(10, 20, 30)) +
  ylab(paste("Number of sex-biased splicing events"))


ggsave(file = "fig2d.pdf", plot = p)

