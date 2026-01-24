#!/usr/bin/env python3
"""
Download SHA order 3 equation data from Tom Fisher's website.

Downloads 30 files from:
https://www.dpmms.cam.ac.uk/~taf1000/g1data/sha3eqns.{A}-{B}

and saves them to data/sha_order3_input/sha3eqns.{A}-{B}.txt
"""

import os
import urllib.request
import ssl

BASE_URL = "https://www.dpmms.cam.ac.uk/~taf1000/g1data/"
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "data", "sha_order3_input")

# Generate the 30 file ranges
# First file is 00000-09999, then 10000-19999, ..., 290000-299999
FILE_RANGES = []
for i in range(30):
    start = i * 10000
    end = start + 9999
    if i == 0:
        # First range has leading zeros
        FILE_RANGES.append(f"{start:05d}-{end:05d}")
    else:
        FILE_RANGES.append(f"{start}-{end}")


def download_files():
    """Download all 30 SHA3 equation files."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Create SSL context that doesn't verify certificates
    # (the site has certificate issues)
    ssl_context = ssl.create_default_context()
    ssl_context.check_hostname = False
    ssl_context.verify_mode = ssl.CERT_NONE

    for file_range in FILE_RANGES:
        filename = f"sha3eqns.{file_range}"
        url = f"{BASE_URL}{filename}"
        output_path = os.path.join(OUTPUT_DIR, f"{filename}.txt")

        if os.path.exists(output_path):
            print(f"Skipping {filename} (already exists)")
            continue

        print(f"Downloading {filename}...")
        try:
            request = urllib.request.Request(url)
            with urllib.request.urlopen(request, context=ssl_context) as response:
                content = response.read()
                with open(output_path, 'wb') as f:
                    f.write(content)
            print(f"  Saved to {output_path}")
        except Exception as e:
            print(f"  Error downloading {filename}: {e}")

    print("\nDownload complete!")


if __name__ == "__main__":
    download_files()
