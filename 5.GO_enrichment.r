# GO enricher
library(clusterProfiler)
input <- read.delim("inputID.txt", row.names=NULL)
df <- enricher(gene=input,
             pvalueCutoff = 0.05,
             pAdjustMethod = "BH",
             TERM2GENE = term2gene2,
             TERM2NAME = term2name,
             universe = universe)
df <- as.data.frame(df)
write.table(df, file="df.txt",  sep="\t")

# Plot network
install.packages("aPEAR")
library(aPEAR)

df$padj = -log10(df$p.adjust)
enrichmentNetwork(
  df,
  colorBy = 'padj',
  nodeSize = 'Count' ,
  drawEllipses = TRUE
)+
  scale_color_gradientn(colours = c("#A5A6CF",'#C9BC9C','#66C6F2'),
                        name = "-log10(p.adjust)")