#Libraries___________________________________________________________________________________________________________________________________________________________________________________________________________________
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)
packageVersion("ComplexHeatmap")

#____________________________________________________________________________________________________________________________________________________________________________________________

# ---------------------------
# Read & prepare data
# ---------------------------
genes <- read_excel("Final_table.xlsx",.name_repair = "minimal")
fpkm<- genes[, grep("FPKM", colnames(genes))]
fpkm$Gene_id<-genes$Gene_id

# Convert wide to long format
df_long <- fpkm %>%
  pivot_longer(
    cols = -Gene_id,
    names_to = "Sample",
    values_to = "FPKM"
  )

# Log10 transform
df_long$logFPKM <- log10(df_long$FPKM + 0.001)

# Density plot____________________________________________________________________________________________________________________________________________________________________________________________
ggplot(df_long, aes(x = logFPKM, color = Sample, fill = Sample)) +
  geom_density(alpha = 0.15, size = 0.5) +
  labs(x = expression(bold(log[10](FPKM))), y = "Density"
  ) +
  theme_bw(base_size = 14) +
  coord_cartesian(ylim = c(0, 0.65)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 0.6, by = 0.2)
  ) +
  scale_color_discrete(
    labels = c("V592-1", "V592-2", "V592-3",
      expression(Delta*bolditalic("VdPT1-1")),
      expression(Delta*bolditalic("VdPT1-2")),
      expression(Delta*bolditalic("VdPT1-3"))
    )
  ) +
  scale_fill_discrete(
    labels = c("V592-1","V592-2","V592-3",
      expression(Delta*bolditalic("VdPT1-1")),
      expression(Delta*bolditalic("VdPT1-2")),
      expression(Delta*bolditalic("VdPT1-3"))
    )
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(hjust = 1, size = 18, face = "bold", color = "black",family = "serif"),
    axis.text.x = element_text(hjust = 1, size = 18, face = "bold", color = "black", family = "serif"),
    axis.title.y = element_text(size = 20, face = "bold", color = "black", family = "serif"),
    axis.title.x = element_text(size = 20, face = "bold", color = "black", family = "serif"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = c(0.95, 0.95),
    legend.justification = c("right", "top"),
    legend.background = element_rect(color = "transparent", fill = "transparent", size = 0.3),
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(face = "bold", family = "serif", size = 16),
    legend.title = element_blank()
  )


ggsave(
  filename = "density_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 6,                             # width in inches
  height = 5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)

# Boxplot____________________________________________________________________________________________________________________________________________________________________________________________

ggplot(df_long, aes(x = Sample, y = logFPKM, fill = Sample, color = Sample)) +
  
  geom_boxplot(outlier.size = 0.5, linewidth = 0.7, alpha = 0.6) +
  
  theme_bw() +
  
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  
  labs(x = "",y = expression(bold(log[10](FPKM)))) +
  theme(
    axis.text.x = element_text(angle = 45, color = "black",hjust = 1, face = "bold", family = "serif", size =18),
    axis.text.y = element_text(color = "black", size =18, family = "serif", face = "bold"),
    axis.title = element_text(color = "black", size =20, family = "serif", face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )+
  scale_x_discrete(labels = c(
    "V592-1_FPKM" = expression(bold("V592-1")),
    "V592-2_FPKM" = expression(bold("V592-2")),
    "V592-3_FPKM" = expression(bold("V592-3")),
    "∆VdPT1-1_FPKM" = expression(Delta*bolditalic("VdPT1-1")),
    "∆VdPT1-2_FPKM" = expression(Delta*bolditalic("VdPT1-2")),
    "∆VdPT1-3_FPKM" = expression(Delta*bolditalic("VdPT1-3"))
  ))
  
ggsave(
  filename = "box_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 6,                             # width in inches
  height = 5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)

#volcano plot____________________________________________________________________________________________________________________________________________________________________________________________

df <- genes[,c(5,9)]
df$log2FoldChange <- as.numeric(as.character(df$log2FoldChange))
df$padj <- as.numeric(as.character(df$padj))
df <- na.omit(df)

df$significant <- ifelse(df$padj < 0.05 & df$log2FoldChange >= 1, "Up",
                         ifelse(df$padj < 0.05 & df$log2FoldChange <= -1, "Down",
                                "Normal"))

ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.7, size = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(x = expression(bold(log[2](FC))), y = expression(bold(-log[10](FDR)))) +
  scale_color_manual(values = c(
    "Down" = "#2C7BB6",
    "Normal" = "#7A7A7A",
    "Up" = "#D7191C"
  )) +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = c(0.9, 0.87),
        axis.text.x = element_text(color = "black",hjust = 1, face = "bold", family = "serif", size =18),
        axis.text.y = element_text(color = "black", size =18, family = "serif", face = "bold"),
        axis.title = element_text(color = "black", size =20, family = "serif", face = "bold"),
        legend.text = element_text(face = "bold", family = "serif", size = 16))

ggsave(
  filename = "volcano_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 6.5,                             # width in inches
  height = 4.5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)

# MA plot____________________________________________________________________________________________________________________________________________________________________________________________________________________

fpkm_cols <- genes[, grep("FPKM", colnames(genes))]
meanFPKM <- rowMeans(fpkm_cols, na.rm = TRUE)

df1 <- data.frame(
  log2FoldChange = as.numeric(genes$log2FoldChange),
  meanFPKM = meanFPKM,
  padj = as.numeric(genes$padj)
)
df1 <- na.omit(df1)
df1$significant <- ifelse(df1$padj < 0.05 & df1$log2FoldChange >= 1, "Up",
                          ifelse(df1$padj < 0.05 & df1$log2FoldChange <= -1, "Down",
                                 "Normal"))

df1$significant <- factor(df1$significant, levels = c("Down", "Normal", "Up"))

# MA plot
ggplot(df1, aes(x = log2(meanFPKM + 1), y = log2FoldChange)) +
  geom_point(aes(color = significant), size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = expression(bold(log[10](FPKM))), y = expression(bold(log[2](FC))))+
  scale_color_manual(values = c("Up" = "#D7191C","Normal" = "#7A7A7A","Down" = "#2C7BB6"
  )) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = c(0.85, 0.87),
        axis.text.x = element_text(color = "black",hjust = 1, face = "bold", family = "serif", size =18),
        axis.text.y = element_text(color = "black", size =18, family = "serif", face = "bold"),
        axis.title = element_text(color = "black", size =20, family = "serif", face = "bold"),
        legend.text = element_text(face = "bold", family = "serif", size = 16))

ggsave(
  filename = "MA_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 6.5,                             # width in inches
  height = 4.5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)


#PCA_________________________________________________________________________________________________________________________________________________________________________________________________________

fpkm <- genes[, grep("FPKM", colnames(genes))]
fpkm <- data.frame(lapply(fpkm, as.numeric))
fpkm_log <- log2(fpkm + 1)
fpkm_log <- fpkm_log[complete.cases(fpkm_log), ]
fpkm_log <- fpkm_log[apply(fpkm_log, 1, sd) > 0, ]
pca <- prcomp(t(fpkm_log), scale. = TRUE)
pca_df <- data.frame(Sample = colnames(fpkm_log), PC1 = pca$x[,1], PC2 = pca$x[,2])
pca_df$Group <- ifelse(grepl("V592", pca_df$Sample), "V592",
                       ifelse(grepl("VdPT1", pca_df$Sample), "ΔVdPT1", NA))
var_exp <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(aes(shape = Group), size = 7) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2, color = NA) +
  theme_bw(base_size = 14) +
  labs(x = paste0("PC1 (", round(var_exp[1], 1), "%)"),
       y = paste0("PC2 (", round(var_exp[2], 1), "%)")) +
  scale_color_discrete(labels = c("V592", expression(Delta*bolditalic("VdPT1")))) +
  scale_fill_discrete(labels = c("V592", expression(Delta*bolditalic("VdPT1")))) +
  scale_shape_discrete(labels = c("V592", expression(Delta*bolditalic("VdPT1")))) +
  guides(fill = "none") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        legend.position = c(0.7, 0.1),
        axis.text.x = element_text(color = "black", hjust = 1, face = "bold", family = "serif", size = 18),
        axis.text.y = element_text(color = "black", size = 18, family = "serif", face = "bold"),
        axis.title = element_text(color = "black", size = 20, family = "serif", face = "bold"),
        legend.text = element_text(face = "bold", family = "serif", size = 18))

ggsave(
  filename = "PCA_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 4.3,                             # width in inches
  height = 5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)



#Hierarchial clustering heatmap____________________________________________________________________________________________________________________________________________________________________________________________

genes$log2FoldChange <- as.numeric(genes$log2FoldChange)
genes$padj <- as.numeric(genes$padj)
deg <- genes[genes$padj < 0.05 & abs(genes$log2FoldChange) >= 1, ]
fpkm <- deg[, grep("FPKM", colnames(deg))]
fpkm <- as.data.frame(lapply(fpkm, as.numeric))
fpkm_log <- log2(fpkm + 1)
fpkm_log <- na.omit(fpkm_log)
fpkm_log <- fpkm_log[apply(fpkm_log, 1, sd) > 0, ]
heatmap_matrix <- t(scale(t(fpkm_log)))
heatmap_matrix <- na.omit(heatmap_matrix)
heatmap_matrix
new_names <- c("V592-2", "V592-1", "ΔVdPT1-3",
               "ΔVdPT1-2", "V592-3", "ΔVdPT1-1")
colnames(heatmap_matrix) <- new_names

png("hierarchial_clustering.png",
    width = 3, height = 8, units = "in", res = 300)

format_labels_fc <- function(column_names) {
  sapply(column_names, function(x) {
    if (grepl("_", x)) {
      gene <- sub("_.*", "", x)
      rest <- sub("^[^_]+", "", x)
      paste0("bold(bolditalic('", gene, "')*'", rest, "')")
    } else {
      paste0("bold('", x, "')")
    }
  })
}

heatmap2 <- Heatmap(
  heatmap_matrix,
  col = colorRampPalette(c("#4575b4", "#f7f7b9", "#d73027"))(1000),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,  # CHANGED TO FALSE
  show_column_names = TRUE,
  column_labels = parse(text = format_labels_fc(colnames(heatmap_matrix))),
  column_names_rot = 45,
  width = ncol(heatmap_matrix) * unit(20, "pt"),
  height = nrow(heatmap_matrix) * unit(0.2, "pt"),
  row_names_gp = gpar(fontfamily = "serif", fontsize = 10, fontface = "bold"),
  column_names_gp = gpar(fontfamily = "serif", fontsize = 12, fontface = "bold"),
  rect_gp = gpar(col = "white", lwd = 0),
  heatmap_legend_param = list(
    title = "FPKM",
    title_gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold"),
    labels_gp = gpar(fontsize = 12, fontfamily = "serif", fontface = "bold"),
    at = seq(-1, 1, 0.5),
    labels = c("-2.00", "-1.00", "0.00", "1.00", "2.00"),
    legend_height = unit(4, "cm"),
    legend_width = unit(6, "cm"),
    title_position = "topcenter",
    border = FALSE
  )
)
draw(heatmap2)

dev.off()



# Correlation matrix____________________________________________________________________________________________________________________________________________________________________________________________


fpkm <- genes[, grep("FPKM", colnames(genes))]
fpkm <- as.data.frame(lapply(fpkm, as.numeric))

fpkm_log <- log2(fpkm + 1)
fpkm_log <- na.omit(fpkm_log)
fpkm_log <- fpkm_log[apply(fpkm_log, 1, sd) > 0, ]
colnames(fpkm_log) <- c("V592-2", "V592-1", "ΔVdPT1-3","ΔVdPT1-2",
                        "V592-3","ΔVdPT1-1" )

cor_matrix <- cor(fpkm_log)


cor_matrix <- cor(fpkm_log, use = "pairwise.complete.obs")

library(pheatmap)

my_labels <- c(
  expression(bold("V592-2")),
  expression(bold("V592-1")),
  expression(Delta * bolditalic("VdPT1-3")),
  expression(Delta * bolditalic("VdPT1-2")),
  expression(bold("V592-3")),
  expression(Delta * bolditalic("VdPT1-1"))

  
)

pheatmap(
  cor_matrix,
  breaks = seq(0.9, 1, by = 0.001),
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 12,
  
)
# Create fine-grained color palette (100+ colors for 0.001 precision)
pheatmap(
  cor_matrix,
  breaks = seq(0.9, 1, by = 0.001),  # Break every 0.001
  display_numbers = TRUE,
  number_format = "%.2f",  # Show 2 decimals
  clustering_distance_rows = "euclidean",
  fontfamily = "serif",
  fontface = "bold",
  fontsize = 12,
  name = " ",
  labels_row = my_labels,
  labels_col = my_labels,
  legend = FALSE)
  


ggsave(
  filename = "correlation_plot.png",  # output file name
  plot = last_plot(),                    # or assign your plot to a variable like p
  width = 6,                             # width in inches
  height = 5,                            # height in inches
  dpi = 300                               # 300 dpi for publication
)


----------------
  deg <- data.frame(
    Regulation = c("Total", "Upregulated", "Downregulated"),
    Count = c(2127, 956, 1171)
  )

deg$Regulation <- factor(deg$Regulation,
                         levels = c("Total", "Upregulated", "Downregulated"))

library(ggplot2)

ggplot(deg, aes(x = Regulation, y = Count, fill = Regulation)) +
  
  geom_bar(stat = "identity", width = 0.6) +
  
  geom_text(
    aes(label = Count),
    vjust = -0.4,
    size = 6,
    family = "serif",
    fontface = "bold"
  ) +
  
  scale_fill_manual(values = c(
    "Total" = "#4D4D4D",
    "Upregulated" = "#D73010",
    "Downregulated" = "#4575B4"
  )) +
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  
  theme_bw() +
  
  labs(
    x = "",
    y = "Number of DEGs"
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      angle = 45,
      color = "black",
      hjust = 1,
      face = "bold",
      family = "serif",
      size = 18
    ),
    
    axis.text.y = element_text(
      color = "black",
      size = 18,
      family = "serif",
      face = "bold"
    ),
    
    axis.title = element_text(
      color = "black",
      size = 20,
      family = "serif",
      face = "bold"
    ),
    
    plot.title = element_text(
      face = "bold",
      family = "serif",
      size = 22,
      hjust = 0.5
    ),
    
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )
