import os
import glob
import hashlib

def check_md5_all(md5_file, fastq_dir):
    # Read MD5 file
    md5db = {}
    with open(md5_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            hash_val, filename = line.split(None, 1)  # split on whitespace
            md5db[filename.strip()] = hash_val.strip()

    mismatches = []

    # Check each FASTQ file
    for f in glob.glob(os.path.join(fastq_dir, "*.fastq*")):
        with open(f, "rb") as fh:
            file_hash = hashlib.md5(fh.read()).hexdigest()

        filename = os.path.basename(f)
        expected_hash = md5db.get(filename)

        if expected_hash != file_hash:
            mismatches.append((filename, expected_hash, file_hash))

    # Report results
    if mismatches:
        print("❌ MD5 mismatches found:")
        for filename, expected, computed in mismatches:
            print(f"{filename} → Expected: {expected}, Computed: {computed}")
    else:
        print("✅ All FASTQ files passed MD5 check")
