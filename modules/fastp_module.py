import os, glob, subprocess
from config import RAW_DIR, TRIMMED_DIR, THREADS

os.makedirs(TRIMMED_DIR, exist_ok=True)

def run_fastp():
    r1_files = sorted(glob.glob(os.path.join(RAW_DIR, "*_1.fastq.gz")))
    for r1 in r1_files:
        r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")
        sample = os.path.basename(r1).replace("_1.fastq.gz", "")
        out_r1 = os.path.join(TRIMMED_DIR, f"{sample}_1.trimmed.fastq.gz")
        out_r2 = os.path.join(TRIMMED_DIR, f"{sample}_2.trimmed.fastq.gz")
        subprocess.run([
            "fastp",
            "-i", r1,
            "-I", r2,
            "-o", out_r1,
            "-O", out_r2,
            "--qualified_quality_phred", "20",
            "--length_required", "30",
            "--thread", str(THREADS)
        ])
        print(f"✅ Finished trimming {sample}")
