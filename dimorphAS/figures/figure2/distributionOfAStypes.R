## Create Figure 2e -- distribution of the five alternative splicing types in tissues with > 10 AS events

# Removal all variables from workspace
rm(list=ls())

library(viridis)
library(ggplot2)
library(magrittr)
library(dplyr)

# set current directory to file location before running!
df <- read.table(file = "figures/figure2/Significant_events.tsv", sep = "\t", quote = NULL, header = T)
totals <- read.table("figures/figure1/Totals_by_tissue_annotated.txt", sep = "\t", quote = NULL, header = T)


idx <- match(df$Tissue, totals$Tissue)
df$Label <- totals$Label[idx]
df$Label <- factor(df$Label, levels = totals$Label)

head(df, 10)

tissues_keep <- subset(totals, Total > 10)

df_plot <- data.frame()

for (i in 1:length(tissues_keep$Label)) {
  df_tissue <- df %>% dplyr::filter(Label %in% tissues_keep$Label[i])
  df_tissue$Perc <- (df_tissue$Counts / sum(df_tissue$Counts)) * 100
  df_plot <- rbind(df_plot, df_tissue)
}  

df_plot$Label <- factor(df_plot$Label, levels = totals$Label)

p<-ggplot(df_plot, 
       aes(x=Label, y=Perc, fill = ASE, by = Tissue)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, angle = 270, hjust = 0.0, vjust = 0.5),
        #axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(face="plain", colour="black", size=12),
        axis.title.y = element_text(face="plain", colour="black", size=12),
        legend.title=element_blank(),
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        legend.text = element_text(face="plain", colour="black", size=12)
  ) +
  #scale_fill_brewer(palette="Dark2") +
  #scale_fill_manual(values = palette_npg) +
  scale_fill_viridis(discrete=TRUE) +
  ylab("Splicing type (%)") +
  xlab("Tissue")



ggsave(file = "fig2e.pdf", plot = p)
