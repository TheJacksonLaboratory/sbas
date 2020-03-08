library(dplyr)
library(ggplot2)

#Parse files to create a data frame with counts

files <- list.files(path = "significant_events/", pattern = "*.txt")
as_types <- c("a3ss", "a5ss", "mxe", "ri", "se")

files_aux <- gsub(pattern = ".txt", replacement = "", x = files)
files_aux <- gsub(pattern = "a3ss$|a5ss$|mxe$|ri$|se$", replacement = "", files_aux)

a3ss_annot <- read.table(file = "fromGTF.A3SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
a5ss_annot <- read.table(file = "fromGTF.A5SS.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
mxe_annot <- read.table(file = "fromGTF.MXE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
ri_annot <- read.table(file = "fromGTF.RI.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)
se_annot <- read.table(file = "fromGTF.SE.txt", sep = "\t", quote = "\"", header = T, stringsAsFactors = F)

gene_as <- data.frame()

for (i in 1:length(files)) {
  
  lines <- readLines(paste0("./significant_events/", files[i]))
  
  if(length(lines) > 1){ #has significant events
    events <- read.table(paste0("./significant_events/", files[i]), sep = "\t", skip = 1) 
    
    if(grepl("a3ss.txt$", files[i])){
      idx <- match(events$V1, a3ss_annot$ID)
      res <- data.frame(Tissue = files_aux[i], ASE = "A3SS", 
                        GeneSymbol = a3ss_annot$geneSymbol[idx],
                        chr = a3ss_annot$chr[idx])
    }
    if(grepl("a5ss.txt$", files[i])){
      idx <- match(events$V1, a5ss_annot$ID)
      res <- data.frame(Tissue = files_aux[i], ASE = "A5SS", 
                        GeneSymbol = a5ss_annot$geneSymbol[idx],
                        chr = a5ss_annot$chr[idx])
    }
    if(grepl("mxe.txt$", files[i])){
      idx <- match(events$V1, mxe_annot$ID)
      res <- data.frame(Tissue = files_aux[i], ASE = "MXE", 
                        GeneSymbol = mxe_annot$geneSymbol[idx],
                        chr = mxe_annot$chr[idx])
    }
    if(grepl("se.txt$", files[i])){
      idx <- match(events$V1, se_annot$ID)
      res <- data.frame(Tissue = files_aux[i], ASE = "SE", 
                        GeneSymbol = se_annot$geneSymbol[idx],
                        chr = se_annot$chr[idx])
    }
    if(grepl("ri.txt$", files[i])){
      idx <- match(events$V1, ri_annot$ID)
      res <- data.frame(Tissue = files_aux[i], ASE = "RI", 
                        GeneSymbol = ri_annot$geneSymbol[idx],
                        chr = ri_annot$chr[idx])
    }
    
    gene_as <- rbind(gene_as, res)
    
  } #if has sig. events
  
} #for all files

# Count most frequent spliced genes
res <- gene_as %>% group_by(GeneSymbol) %>% count(GeneSymbol) %>% arrange(desc(n)) %>% as.data.frame()
res$GeneSymbol <- factor(res$GeneSymbol, levels = res$GeneSymbol)
length(res$GeneSymbol)

#Add number of tissues
nTissues <- rep(NA, length(res))
for (i in 1:nrow(res)) {
  df_gene <- gene_as %>% filter(GeneSymbol == res$GeneSymbol[i])
  nTissues[i] <- length(unique(df_gene$Tissue))
}
res$Tissues <- nTissues

# Genes with more than 10 splicing events
ggplot(res[res$n > 10, ], aes(x = GeneSymbol, y = n)) +
  geom_point(size = 4, aes(fill = Tissues, color = Tissues)) +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle = 270, hjust = 0.0, vjust = 0.5),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(face="plain", colour="black", 
                                    size=14),
        axis.title.y = element_text(face="plain", colour="black", 
                                    size=14),
        legend.title=element_blank(),
        legend.text = element_text(face="plain", colour="black", 
                                   size=12)) +
  scale_fill_viridis_c(aesthetics = c("colour", "fill"),
                       option = "plasma",
                       limits = c(1, 30), breaks = c(10, 20, 30)) +
  ylab(paste("Number of sex-biased splicing events")) +
  xlab("Genes")



#Pie chart

#Tissue specific 1 tissue
#Tissue group 2-5 tissues
#Recurrent > 5 tissues

#Pie chart - Number of patients with 0, 1, >1 events
counts <- c(res %>% filter(Tissues == 1) %>% count() %>% as.numeric(),
            res %>% filter(Tissues > 1 & Tissues < 5) %>% count() %>% as.numeric(),
            res %>% filter(Tissues > 5) %>% count() %>% as.numeric())

# Define some colors ideal for black & white print
colors <- c("white","grey70","grey90","grey50","black")
colors <- c("grey90","grey50", "black")

# Calculate the percentage for each category

counts_labels <- round(counts/sum(counts) * 100, 1)
counts_labels <- paste(counts_labels, "%", sep="")
counts_labels <- paste(counts_labels, counts, sep =" ")


# Create a pie chart with defined heading and custom colors
# and labels
pie(counts, main="", col=colors, labels=counts_labels,
    cex=1.2)

# Create a legend at the right   
legend(1.5, 0.5, c("1 tissue","2-5 tissues","> 5 tissues"), cex=1.2, 
       fill=colors)


# Number of splicing events normalized per chromosome
library(rtracklayer)
gencode <- import.gff("/Users/dveiga/tools/Gencode/gencode.v25.annotation.gtf")

exons <- gencode[ gencode$type == "exon", ]
exons <- as.data.frame(exons)

all_chr <- as.character(unique(gene_as$chr))
chr_counts <- rep(0, length(all_chr))

for (i in 1:length(all_chr)) {
  chr_counts[i] <- nrow(exons[exons$seqnames == all_chr[i], ])
}

exon_counts <- data.frame(chr = all_chr, counts = chr_counts)

# Count most frequent spliced chromosomes
res <- gene_as %>% group_by(chr) %>% count(chr) %>% arrange(desc(n)) %>% as.data.frame()
res$chr <- factor(res$chr, levels = res$chr)

idx <- match(res$chr, exon_counts$chr)

res$ExonCounts <- exon_counts$counts[idx]

res$Index <- (res$n / res$ExonCounts) * 1000

res_sorted <- res %>% arrange(desc(Index))
res_sorted$chr <- factor(res_sorted$chr, levels = res_sorted$chr)

ggplot(res_sorted, aes(x = chr, y = Index, size = n)) +
  geom_point(color = "red") +
  theme_bw() +
  theme(axis.text.x = element_text(size=10, angle = 270, hjust = 0.0, vjust = 0.5),
        axis.text.y = element_text(size=16),
        axis.title.x = element_text(face="plain", colour="black", 
                                    size=14),
        axis.title.y = element_text(face="plain", colour="black", 
                                    size=14),
        legend.title=element_blank(),
        legend.text = element_text(face="plain", colour="black", 
                                   size=12)) +
  scale_fill_viridis_c(aesthetics = c("colour", "fill"),
                       option = "plasma",
                       limits = c(1, 650)) +
  ylab(paste("Normalized Number of sex-biased ASE")) +
  xlab("Chromosomes") +
  guides(size = guide_legend(title = "Number of ASE"))

write.table(res_sorted, file = "SplicingIndex_chr.tsv", sep = "\t", quote = F,
            row.names = F)
