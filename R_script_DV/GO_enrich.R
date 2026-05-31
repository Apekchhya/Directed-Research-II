library(topGO)
packageVersion("topGO")
library(readxl)

genes <- read_excel("Final_table.xlsx", .name_repair = "minimal")  # load data

gene_ids <- as.character(genes$Gene_id)  # extract gene IDs

genes$log2FoldChange <- as.numeric(genes$log2FoldChange)  # ensure numeric
genes$padj <- as.numeric(genes$padj)  # ensure numeric
genes$stat <- as.numeric(genes$stat)  # ensure numeric

gene_scores <- genes$stat  # gene scoring vector
names(genes_scores) <- genes$Gene_id  # assign names
gene_scores <- gene_scores[!is.na(gene_scores)]  # remove NA scores

deg_genes <- genes$Gene_id[genes$padj < 0.01 & genes$log2FoldChange <=-1]  # DEG selection

deg_criteria <- expression(padj < 0.01 & log2FoldChange <= -1)  # store criteria

gene2GO <- strsplit(as.character(genes$GO_annotation), ";")  # split GO terms
names(gene2GO) <- genes$Gene_id  # assign gene names
gene2GO <- gene2GO[!is.na(names(gene2GO)) & lengths(gene2GO) > 0] # remove invalid entries

geneList <- as.integer(gene_ids %in% deg_genes)  # 0/1 DEG vector
names(geneList) <- gene_ids  # assign gene names
geneList <- geneList[!is.na(names(geneList))]  # safety cleanup

geneSel <- function(allScore) { return(allScore == 1) }  # define DEG selector

GO_data_BP <- new("topGOdata",
               ontology = "BP",  # Biological Process
               allGenes = geneList,  # named 0/1 vector
               geneSel = geneSel,  # IMPORTANT FIX
               annot = annFUN.gene2GO,  # GO annotation
               gene2GO = gene2GO,
               nodeSize = 10)  # GO mapping
resultFisher_BP <- runTest(GO_data_BP,
                        algorithm = "weight01",  # classic topology
                        statistic = "fisher")  # Fisher test

GO_results_BP <- GenTable(GO_data_BP,
                       classicFisher = resultFisher_BP,
                       orderBy = "classicFisher",
                       topNodes = 20)  # top enriched GO terms
GO_results_BP$type <-c('BP')
GO_results_BP


GO_data_MF <- new("topGOdata",
                  ontology = "MF",  # Biological Process
                  allGenes = geneList,  # named 0/1 vector
                  geneSel = geneSel,  # IMPORTANT FIX
                  annot = annFUN.gene2GO,  # GO annotation
                  gene2GO = gene2GO,
                  nodeSize = 10)  # GO mapping
resultFisher_MF <- runTest(GO_data_MF,
                        algorithm = "weight01",  # classic topology
                        statistic = "fisher")  # Fisher test


GO_results_MF <- GenTable(GO_data_MF,
                          classicFisher = resultFisher_MF,
                          orderBy = "classicFisher",
                          topNodes = 20)  # top enriched GO terms
GO_results_MF$type <-c('MF')
GO_results_MF


GO_data_CC <- new("topGOdata",
                  ontology = "CC",
                  allGenes = geneList,
                  geneSel = geneSel,
                  annot = annFUN.gene2GO,
                  gene2GO = gene2GO,
                  nodeSize = 10)

resultFisher_CC <- runTest(GO_data_CC,
                           algorithm = "weight01",
                           statistic = "fisher")

GO_results_CC <- GenTable(GO_data_CC,
                          classicFisher = resultFisher_CC,
                          orderBy = "classicFisher",
                          topNodes = 20)


GO_results_CC$type <- "CC"

GO_results_CC

library(openxlsx)

# 1. Add ontology labels (important for tracking)
GO_results_BP$type <- "BP"
GO_results_CC$type <- "CC"
GO_results_MF$type <- "MF"

# 2. Combine all data frames
GO_all <- rbind(
  GO_results_BP,
  GO_results_CC,
  GO_results_MF
)

# 3. (optional) reorder columns nicely
GO_all <- GO_all[, c("type", names(GO_all)[names(GO_all) != "type"])]

# 4. Export to Excel
write.xlsx(GO_all,
           file = "GO_enrichment_results_down.xlsx",
           rowNames = FALSE)

# 5. check
GO_all
