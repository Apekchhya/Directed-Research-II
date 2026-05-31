import os, subprocess, re
from config import RAW_DIR

def fetch_fastq_from_sh(script_path="/Users/apekshya/Documents/RNA_SEQ/data/data_download.sh"):
    """
    Reads a .sh file, extracts FASTQ URLs, checks which files are missing,
    and downloads them into RAW_DIR using wget.
    """

    os.makedirs(RAW_DIR, exist_ok=True)
    download_dir = RAW_DIR

    # 1. Read URLs from the script
    with open(script_path, "r") as f:
        lines = f.readlines()

    urls = []
    pattern = r"wget\s+(?:-nc\s+)?(ftp\S+\.fastq\.gz)"

    for line in lines:
        match = re.search(pattern, line)
        if match:
            urls.append(match.group(1))

    print(f"Found {len(urls)} URLs in script.")

    # 2. Check which files are missing
    missing_urls = []

    for url in urls:
        filename = url.split("/")[-1]
        filepath = os.path.join(download_dir, filename)

        if not os.path.exists(filepath):
            missing_urls.append(url)

    print(f"Missing files: {len(missing_urls)} / {len(urls)}")

    if len(missing_urls) == 0:
        print("All files already downloaded.")
        return  # exit function

    # 3. Download missing files
    print("Downloading missing files...\n")

    for url in missing_urls:
        try:
            subprocess.run(
                ["wget", "-nc", "-P", download_dir, url],
                check=True
            )
            print(f"Downloaded: {url.split('/')[-1]}")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {url}: {e}")

    print("Done!")
