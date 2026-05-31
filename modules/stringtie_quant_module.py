import os
import glob
import subprocess
from config import HISAT2_OUT, MERGE_OUT, EXPR_OUT, THREADS, PREPDE, ANNOTATION_GTF

os.makedirs(EXPR_OUT, exist_ok=True)

def run_quantification_and_prepde():
    merged_gtf = ANNOTATION_GTF

    if not os.path.exists(merged_gtf):
        raise FileNotFoundError(f"Merged GTF not found: {merged_gtf}")

    bam_files = sorted(glob.glob(f"{HISAT2_OUT}/*.sorted.bam"))
    if not bam_files:
        raise FileNotFoundError(f"No sorted BAM files found in {HISAT2_OUT}")

    print(f"Found {len(bam_files)} BAM files. Running StringTie quantification...")

    for bam in bam_files:
        sample = os.path.basename(bam).replace(".sorted.bam", "")
        sample_dir = os.path.join(EXPR_OUT, sample)
        os.makedirs(sample_dir, exist_ok=True)

        out_gtf = os.path.join(sample_dir, f"{sample}.gtf")

        if os.path.exists(out_gtf):
            print(f"⏭️  Skipping {sample}, GTF already exists")
            continue

        print(f"🔄 Running StringTie for {sample}...")
        subprocess.run([
            "stringtie", bam,
            "-e", "-B",
            "-G", merged_gtf,
            "-o", out_gtf,
            "-p", str(THREADS)
        ], check=True)
        print(f"✅ Quantification done: {sample}")

    # Run prepDE.py3 to generate count matrices
    gene_matrix = os.path.join(EXPR_OUT, "gene_count_matrix.csv")
    transcript_matrix = os.path.join(EXPR_OUT, "transcript_count_matrix.csv")

    if os.path.exists(gene_matrix):
        print("⏭️  Count matrices already exist, skipping prepDE")
        return

    print("🔄 Running prepDE.py3 to generate count matrices...")
    subprocess.run([
        "python3", PREPDE,
        "-i", EXPR_OUT,
        "-g", gene_matrix,
        "-t", transcript_matrix,
        "-l", "150" 
    ], check=True)

    print(f"✅ Count matrices saved in: {EXPR_OUT}")
    print(f"   → gene_count_matrix.csv      (use this for DESeq2)")
    print(f"   → transcript_count_matrix.csv (ignore for DESeq2)")