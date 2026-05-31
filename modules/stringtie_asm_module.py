import os, glob, subprocess
from config import HISAT2_OUT, STRINGTIE_OUT, ANNOTATION_GTF, THREADS

os.makedirs(STRINGTIE_OUT, exist_ok=True)

def run_stringtie():
    bam_files = sorted(glob.glob(f"{HISAT2_OUT}/*.sorted.bam"))

    if not bam_files:
        print("⚠️  No BAM files found in HISAT2_OUT")
        return

    for bam in bam_files:
        sample = os.path.basename(bam).replace(".sorted.bam", "")
        out_gtf = f"{STRINGTIE_OUT}/{sample}.gtf"

        # Skip if already processed
        if os.path.exists(out_gtf):
            print(f"⏭️  Skipping {sample}, GTF already exists")
            continue

        print(f"🔄 Running StringTie for {sample}...")

        subprocess.run([
            "stringtie", bam,
            "-G", ANNOTATION_GTF,
            "-o", out_gtf,
            "-p", str(THREADS),
            "-l", "NewGene",                                        # novel gene label     
            "-C", f"{STRINGTIE_OUT}/{sample}_coverage.gtf",         # coverage GTF
        ], check=True)

        print(f"✅ StringTie assembly done: {sample}")
