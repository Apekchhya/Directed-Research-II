import subprocess
import glob
import sys

# Path to your RNA_SEQ folder
folder_path = "/Users/apekshya/Documents/RNA_SEQ/data/fastq"

# Get all .fa, .fasta, .fq, .fastq files including gzipped versions
files = glob.glob(folder_path + "/*.f[aq]*") + glob.glob(folder_path + "/*.f[aq]*.gz")

if not files:
    print("No FASTA/FASTQ files found in the folder.")
    sys.exit()

# Function to run seqkit stats
def seqkit_stats(file_list):
    try:
        # Build the command: seqkit stats file1 file2 ... -a
        cmd = ["seqkit", "stats"] + file_list + ["-a"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(result.stdout)
    except subprocess.CalledProcessError as e:
        print("Error running seqkit:", e, file=sys.stderr)
        print(e.stdout)
        print(e.stderr)

# Run stats
seqkit_stats(files)