import os, subprocess
from config import MERGE_OUT, GFFCOMPARE_OUT, ANNOTATION_GTF

os.makedirs(GFFCOMPARE_OUT, exist_ok=True)

def run_gffcompare():
    merged_gtf  = f"{MERGE_OUT}/stringtie_merged.gtf"
    output_prefix = f"{GFFCOMPARE_OUT}/gffcompare"

    if not os.path.exists(merged_gtf):
        raise RuntimeError("❌ Merged GTF not found! Run merge step first.")

    if os.path.exists(f"{output_prefix}.tracking"):
        print("⏭️  gffcompare output already exists, skipping")
        return

    print("🔄 Running gffcompare...")

    subprocess.run([
        "gffcompare",
        "-r", ANNOTATION_GTF,      # reference annotation
        "-o", output_prefix,       # output prefix
        "-G",                      # consider only reference transcripts overlapping input
        merged_gtf
    ], check=True)

    print("✅ gffcompare completed successfully")
    print(f"\n📂 Output files:")
    print(f"   {output_prefix}.annotated.gtf  → transcripts with class codes")
    print(f"   {output_prefix}.tmap           → transcript-level classification")
    print(f"   {output_prefix}.refmap         → reference transcript mapping")
    print(f"   {output_prefix}.stats          → summary statistics")
