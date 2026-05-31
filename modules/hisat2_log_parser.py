import re
import pandas as pd
from pathlib import Path
from config import HISAT2_OUT


def generate_hisat2_mapping_table(output_csv):
    """
    Parse HISAT2 log files and generate a mapping summary table.

    Calculations:
    Mapped Reads        = (line4 + line5 + line8) * 2 + line13 + line14
    Uniq Mapped Reads   = (line4 + line8) * 2 + line13
    Multiple Map Reads  = (line5 * 2) + line14
    Total Reads         = nums[0] * 2
    Percentages added in parentheses
    """

    records = []

    # Loop through all HISAT2 log files
    for log_file in Path(HISAT2_OUT).glob("*.err"):
        with open(log_file) as f:
            lines = f.readlines()

        # Extract leading numbers
        nums = []
        for line in lines:
            m = re.match(r"\s*(\d+)", line)
            nums.append(int(m.group(1)) if m else 0)

        # 0-based indexing
        total_reads = nums[0] * 2
        line4 = nums[3]    # concordantly exactly 1 time
        line5 = nums[4]    # concordantly >1 times
        line8 = nums[7]    # discordantly 1 time
        line13 = nums[12]  # mates aligned exactly 1 time
        line14 = nums[13]  # mates aligned >1 times

        # Calculations
        mapped_reads_val = (line4 + line5 + line8) * 2 + line13 + line14
        uniq_mapped_val = (line4 + line8) * 2 + line13
        multi_mapped_val = (line5 * 2) + line14

        # Percentages
        mapped_pct = round((mapped_reads_val / total_reads) * 100, 2)
        uniq_pct = round((uniq_mapped_val / total_reads) * 100, 2)
        multi_pct = round((multi_mapped_val / total_reads) * 100, 2)

        # Format numbers with percentage
        mapped_reads = f"{mapped_reads_val} ({mapped_pct}%)"
        uniq_mapped = f"{uniq_mapped_val} ({uniq_pct}%)"
        multi_mapped = f"{multi_mapped_val} ({multi_pct}%)"

        records.append({
            "Sample": log_file.stem,
            "Total Reads": total_reads,
            "Mapped Reads": mapped_reads,
            "Uniq Mapped Reads": uniq_mapped,
            "Multiple Map Reads": multi_mapped
        })

    # Save to CSV
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"✅ HISAT2 mapping summary saved to: {output_csv}")

    return df
