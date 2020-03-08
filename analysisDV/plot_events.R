library(dplyr)
library(ggplot2)
library(scales)
library(viridis)

totals <- read.table("Totals_by_tissue_annotated.txt", sep = "\t", quote = NULL, header = T)


totals_s <- totals %>% arrange(desc(Total))

totals$Label <- factor(totals$Label, levels = totals_s$Label)

ggplot(totals, aes(x = Label, y = Total, size = Total)) +
  geom_point(color = "red") +
  theme_bw() +
  # scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_continuous(limits = c(0, 140e3), breaks=c(1:6 * 25000)) +
  #coord_trans(y = "log10") +
  scale_y_continuous(trans=log_trans(), breaks=c(1,10,100, 1000, 5000, 10000)) +
  theme(axis.text.x = element_text(size=10, angle = 270, hjust = 0.0, vjust = 0.5),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(face="plain", colour="black", 
                                    size=14),
        axis.title.y = element_text(face="plain", colour="black", 
                                    size=14),
        legend.title=element_blank(),
        legend.text = element_text(face="plain", colour="black", 
                                   size=12)) +
  ylab(paste("Number of sex-biased splicing events")) +
  xlab("Tissue") + 
  guides(size=FALSE)


df <- read.table(file = "Significant_events.tsv", sep = "\t", quote = NULL, header = T)

idx <- match(df$Tissue, totals$Tissue)
df$Label <- totals$Label[idx]
df$Label <- factor(df$Label, levels = totals_s$Label)

head(df, 10)

tissues_keep <- subset(totals, Total > 10)

df_plot <- data.frame()

for (i in 1:length(tissues_keep$Label)) {
  
  df_tissue <- df %>% filter(Label %in% tissues_keep$Label[i])
  df_tissue$Perc <- (df_tissue$Counts / sum(df_tissue$Counts)) * 100
  df_plot <- rbind(df_plot, df_tissue)
    
}  

df_plot$Label <- factor(df_plot$Label, levels = totals_s$Label)

ggplot(df_plot, 
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
