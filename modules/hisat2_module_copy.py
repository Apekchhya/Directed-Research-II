import os, glob, subprocess
from config import RAW_DIR, HISAT2_INDEX_PREFIX, HISAT2_OUT, REFERENCE_FASTA, THREADS

os.makedirs(HISAT2_OUT, exist_ok=True)

def hisat2_index_exists(prefix):
    return all(os.path.exists(f"{prefix}.{i}.ht2") for i in range(1, 9))

def run_hisat2():
    if not hisat2_index_exists(HISAT2_INDEX_PREFIX):
        print("🔨 Building HISAT2 index...")
        subprocess.run(["hisat2-build", REFERENCE_FASTA, HISAT2_INDEX_PREFIX], check=True)
        print("✅ Index built!")

    r1_files = sorted(glob.glob(os.path.join(RAW_DIR, "*_1.fastq.gz")))

    if not r1_files:
        print("⚠️ No FASTQ files found in RAW_DIR")
        return

    for r1 in r1_files:
        r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")
        sample = os.path.basename(r1).replace("_1.fastq.gz", "")
        bam_file = os.path.join(HISAT2_OUT, f"{sample}.sorted.bam")
        sam_file = os.path.join(HISAT2_OUT, f"{sample}.sam")

        # Skip if already done
        if os.path.exists(bam_file):
            print(f"⏭️  Skipping {sample}, BAM already exists")
            continue

        # Check R2 exists
        if not os.path.exists(r2):
            print(f"⚠️  R2 not found for {sample}, skipping")
            continue

        print(f"🔄 Aligning {sample}...")

        # Align with HISAT2
        with open(os.path.join(HISAT2_OUT, f"{sample}.log"), "w") as log, \
            open(os.path.join(HISAT2_OUT, f"{sample}.err"), "w") as err:
            subprocess.run([
                "hisat2", "-p6",
                "--dta",
                "-x", HISAT2_INDEX_PREFIX,
                "-1", r1,
                "-2", r2,
                "-S", sam_file
            ], stdout=log, stderr=err, check=True)

        # SAM → sorted BAM → index
        subprocess.run(
            f"samtools view -b {sam_file} | samtools sort -@ {THREADS} -o {bam_file}",
            shell=True, check=True
        )
        subprocess.run(["samtools", "index", bam_file], check=True)

        # Clean up SAM
        os.remove(sam_file)
        print(f"✅ Finished HISAT2 alignment for {sample}")