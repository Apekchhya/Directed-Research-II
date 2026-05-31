from modules.stringtie_asm_module import run_stringtie
from modules.stringtie_merge_module import run_merge
from modules.gff_compare_module import run_gffcompare
from modules.stringtie_quant_module import run_quantification_and_prepde
from modules.hisat2_module_copy import run_hisat2
from modules.data_integrity_check_module import check_md5_all
from config import EXPR_OUT, MD5_FILE, RAW_DIR, STRINGTIE_OUT
from modules.fetchdata_module import fetch_fastq_from_sh
from modules.hisat2_log_parser import generate_hisat2_mapping_table
from config import HISAT2_OUT
from modules.deseq2_runner import run_deseq2
output_csv = f"{HISAT2_OUT}/hisat2_mapping_summary.csv"




#fetch_fastq_from_sh()
# check_md5_all(MD5_FILE, RAW_DIR)
# run_hisat2()
# generate_hisat2_mapping_table(output_csv)
#run_featurecounts()

# run_stringtie()
# run_merge()
# run_gffcompare()
# run_quantification_and_prepde()
run_deseq2()


# extract_fpkm_matrix(
#     stringtie_dir=EXPR_OUT,
#     output_csv=f"{EXPR_OUT}/fpkm_matrix.csv"
# )