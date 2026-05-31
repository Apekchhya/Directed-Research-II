# Directed-Research-II---Transcriptomics
The repository contains pipeline for RNA seq data analysis aim to reproduce the output of paper entitled "VdPAT1 encoding a pantothenate transporter protein is required for fungal growth, mycelial penetration and pathogenicity of Verticillium dahliae." The analysis work flow is

Raw RNA-seq Data (FASTQ files)
            │
            ▼
Quality Assessment and Data Inspection
(seqkit v2.13.0)
            │
            ▼
Genome Alignment
(HISAT2 v2.2.2)
            │
            ▼
Transcript and Gene Quantification
(StringTie v3.0.3)
(Counts, FPKM and TPM)
            │
            ▼
Differential Expression Analysis
(DEGs Identification)
            │
            ├─────────────┐
            ▼             ▼
GO Enrichment      KEGG Enrichment
(topGO v2.62.0)   (clusterProfiler v4.18.4)
            │             │
            └──────┬──────┘
                   ▼
Biological Interpretation
            │
            ▼
Visualization
(ggplot2 v4.0.3,
ComplexHeatmap v2.26.1)
            │
            ▼
Transcription Factor Analysis
(CIS-BP Database v3.00)
