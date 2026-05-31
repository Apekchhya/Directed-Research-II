import os, glob, subprocess
from config import TRIMMED_DIR, HISAT2_INDEX_PREFIX, HISAT2_OUT,REFERENCE_FASTA, THREADS

os.makedirs(HISAT2_OUT, exist_ok=True)

def hisat2_index_exists(prefix):
    return all(os.path.exists(f"{prefix}.{i}.ht2") for i in range(1, 9))

def run_hisat2():
    if not hisat2_index_exists(HISAT2_INDEX_PREFIX):
        print("🔨 Building HISAT2 index ...")
        subprocess.run(["hisat2-build", REFERENCE_FASTA, HISAT2_INDEX_PREFIX], check=True)
        print("✅ Index built!")
    
    r1_files = sorted(glob.glob(os.path.join(TRIMMED_DIR, "*_1.trimmed.fastq.gz")))
    for r1 in r1_files:
        r2 = r1.replace("_1.trimmed.fastq.gz", "_2.trimmed.fastq.gz")
        sample = os.path.basename(r1).replace("_1.trimmed.fastq.gz", "")
        sam_file = os.path.join(HISAT2_OUT, f"{sample}.sam")
        bam_file = os.path.join(HISAT2_OUT, f"{sample}.sorted.bam")
        # Align
        with open(os.path.join(HISAT2_OUT, f"{sample}.log"), "w") as log:
            subprocess.run([
                "hisat2", "-p", str(THREADS),
                "-x", HISAT2_INDEX_PREFIX,
                "-1", r1,
                "-2", r2,
                "-S", sam_file
            ], stdout=log, stderr=log, check=True)
        # SAM → BAM
        subprocess.run(f"samtools view -b {sam_file} | samtools sort -@ {THREADS} -o {bam_file}", shell=True, check=True)
        subprocess.run(["samtools", "index", bam_file], check=True)
        os.remove(sam_file)
        print(f"✅ Finished HISAT2 alignment for {sample}")
