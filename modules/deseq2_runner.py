import os
import subprocess
from config import EXPR_OUT, DESEQ2_OUT, METADATA_CSV

os.makedirs(DESEQ2_OUT, exist_ok=True)

def run_deseq2():
    gene_matrix = os.path.join(EXPR_OUT, "ref_count.csv")
    r_script    = os.path.join(os.path.dirname(__file__), "deseq2_analysis.R")

    if not os.path.exists(gene_matrix):
        raise FileNotFoundError(f"Count matrix not found: {gene_matrix}\nRun quantification step first.")

    if not os.path.exists(METADATA_CSV):
        raise FileNotFoundError(f"Metadata CSV not found: {METADATA_CSV}")

    results_file = os.path.join(DESEQ2_OUT, "deseq2_genes_results.csv")
    if os.path.exists(results_file):
        print("⏭️  DESeq2 results already exist, skipping")
        return

    print("🔄 Running DESeq2 analysis...")
    subprocess.run([
        "Rscript", r_script,
        gene_matrix,
        METADATA_CSV,
        DESEQ2_OUT
    ], check=True)

    print(f"✅ DESeq2 completed!")
    print(f"   → {DESEQ2_OUT}/deseq2_genes_results.csv")
    print(f"   → {DESEQ2_OUT}/upregulated_genes.csv")
    print(f"   → {DESEQ2_OUT}/downregulated_genes.csv")