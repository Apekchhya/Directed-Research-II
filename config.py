BASE_DIR = "/Users/apekshya/Documents/RNA_SEQ"


#Path configurations
MD5_FILE = f"{BASE_DIR}/data/md5_checksums.txt"
HISAT2_INDEX_PREFIX = f"{BASE_DIR}/data/reference/vd_index"
REFERENCE_FASTA = f"{BASE_DIR}/data/reference/Verticillium_dahliae.ASM15067v2.dna.toplevel.fa"
ANNOTATION_GTF= f"{BASE_DIR}/data/reference/annotation.gtf"
PREPDE = f"{BASE_DIR}/scripts/prepDE.py3"
METADATA_CSV = f"{BASE_DIR}/data/metadata.csv"
MERGE_LABEL = "Newgene"
#output directories
RAW_DIR = f"{BASE_DIR}/data/fastq"
TRIMMED_DIR = f"{BASE_DIR}/results/trimmed_data"
HISAT2_OUT = f"{BASE_DIR}/results/hisat2"
STRINGTIE_OUT = f"{BASE_DIR}/results/stringtie"
MERGE_OUT = f"{BASE_DIR}/results/merge"
GFFCOMPARE_OUT = f"{BASE_DIR}/results/gffcompare"
EXPR_OUT = f"{BASE_DIR}/results/expression"
DESEQ2_OUT = f"{BASE_DIR}/results/deseq2"
THREADS = 4
