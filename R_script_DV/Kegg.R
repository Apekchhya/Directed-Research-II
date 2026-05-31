
# Load library
library(readxl)
library(clusterProfiler)
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)
?clusterProfiler
genes <- read_excel("Final_table.xlsx",.name_repair = "minimal")
tail(genes,5)
universe<- genes$Gene_symbol
universe
genes$log2FoldChange <- as.numeric(genes$log2FoldChange)
genes$padj <- as.numeric(genes$padj)

down_genes<- subset(genes, padj < 0.01 & log2FoldChange <= -1)
up_genes<-subset(genes, padj < 0.01 & log2FoldChange >= 1)

down_genes<-down_genes$Gene_symbol
up_genes<-up_genes$Gene_symbol


kegg_up <- enrichKEGG(
  gene         = up_genes,
  organism     = 'vda',      # Verify this is the code for Verticillium dahliae
  keyType      = 'kegg',     # Use 'kegg' for Locus Tags like VDAG_...
  pvalueCutoff = 5,
  qvalueCutoff = 5
)

kegg_df_up <- as.data.frame(kegg_up)
kegg_df_up
write.xlsx(kegg_df_up, "KEGG_enrichment_up.xlsx", rowNames = FALSE)

kegg_down <- enrichKEGG(
  gene         = down_genes,
  organism     = 'vda',      # Verify this is the code for Verticillium dahliae
  keyType      = 'kegg',     # Use 'kegg' for Locus Tags like VDAG_...
  pvalueCutoff = 5,
  qvalueCutoff = 5
)

kegg_df_down <- as.data.frame(kegg_down)
write.xlsx(kegg_df_down, "KEGG_enrichment_down.xlsx", rowNames= FALSE)
