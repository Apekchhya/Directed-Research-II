import os
import glob
import subprocess
from config import STRINGTIE_OUT, MERGE_OUT, ANNOTATION_GTF, MERGE_LABEL


def run_merge():
    # Create output directory
    os.makedirs(MERGE_OUT, exist_ok=True)

    merged_gtf = os.path.join(MERGE_OUT, "stringtie_merged.gtf")
    mergelist  = os.path.join(MERGE_OUT, "mergelist.txt")

    # Skip if already done
    if os.path.exists(merged_gtf):
        print("⏭️  Merged GTF already exists, skipping")
        return

    # Collect all GTF files
    gtf_files = [
        f for f in glob.glob(f"{STRINGTIE_OUT}/**/*.gtf", recursive=True)
        if not f.endswith("_coverage.gtf")
    ]

    if not gtf_files:
        raise RuntimeError("❌ No GTF files found for merging!")

    print(f"🔎 Found {len(gtf_files)} GTF files for merging")

    # Write mergelist file
    with open(mergelist, "w") as f:
        f.write("\n".join(sorted(gtf_files)))

    print("📝 mergelist.txt created")

    # Run StringTie merge (NO filters)
    print("🔄 Running StringTie merge...")

    subprocess.run([
        "stringtie",
        "--merge",
        "-l", MERGE_LABEL,
        "-G", ANNOTATION_GTF,
        "-o", merged_gtf,
        "-F", "0.01",  # No filtering
        "-T", "0.0",
        "gap-avg", "0.0",
        mergelist
    ], check=True)

    print("✅ StringTie merge completed")
    print(f"📂 Output file: {merged_gtf}")