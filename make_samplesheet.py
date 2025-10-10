#!/usr/bin/env -S uv run --script
#./
# /// script
# requires-python = ">=3.12"
# dependencies = [
# ]
# ///
"""
Make a samplesheet CSV from an Illumina-style FASTQ directory.

Usage:
  ./make_samplesheet.py /home/marcus/projects/mot23/DEMUX

Notes:
- Pairs files matching *_R1_*.fastq.gz with *_R2_* mate.
- Sample name is trimmed before the barcode chunk (e.g. before `_09Z_`).
- Skips files containing '.first' by default (e.g. '.first400k.fastq.gz').
"""

import argparse
import csv
import re
import sys
from pathlib import Path

# Pattern to trim sample name *before* the index/barcode chunk like _09Z_ or _10A_
TRIM_BEFORE_BARCODE = re.compile(r"_[0-9]{2}[A-Z]?_")

def derive_sample_name(filename: str) -> str:
    """Return the sample name before the barcode chunk, else before lane/S/R."""
    import re
    # Look for _S followed by digits (Illumina sample number pattern)
    m = re.search(r'_S\d+', filename)
    if m:
        return filename[:m.start()]
    # Fallback to _L pattern
    idx = filename.find("_L")
    if idx != -1:
        return filename[:idx]
    # As a last resort, strip extension and return
    return filename.replace(".fastq.gz", "").replace(".fastq", "")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Directory containing FASTQ files (default: current directory)",
    )
    ap.add_argument(
        "--include-first",
        action="store_true",
        help="Include files containing '.first' (e.g. '.first400k.fastq.gz').",
    )
    ap.add_argument(
        "-o", "--output",
        default=None,
        help="Output CSV filename (default: samplesheet.csv in the input directory).",
    )
    args = ap.parse_args()

    d = Path(args.directory).expanduser().resolve()
    if not d.is_dir():
        print(f"ERROR: {d} is not a directory", file=sys.stderr)
        sys.exit(1)

    output_csv = Path(args.output) if args.output else d / "samplesheet.csv"

    # Collect FASTQs
    fastqs = sorted(list(d.glob("*.fastq.gz")) + list(d.glob("*.fastq")))
    if not args.include_first:
        fastqs = [p for p in fastqs if ".first" not in p.name]

    # Map R1 -> R2 partners
    r1s = [p for p in fastqs if "_R1_" in p.name or "_R1." in p.name]
    fileset = {p.name for p in fastqs}

    rows = []
    missing_r2 = []
    for r1 in r1s:
        if "_R1_" in r1.name:
            r2_name = r1.name.replace("_R1_", "_R2_")
        else:  # "_R1." pattern
            r2_name = r1.name.replace("_R1.", "_R2.")
        if r2_name not in fileset:
            missing_r2.append(r1.name)
            continue
        sample = derive_sample_name(r1.stem)  # stem removes .gz then .fastq
        rows.append([
            sample,
            str(r1),
            str(d / r2_name),
        ])

    # Write CSV
    with output_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample", "fastq_1", "fastq_2"])
        # sort by sample name for consistency
        for row in sorted(rows, key=lambda r: r[0]):
            w.writerow(row)

    print(f"Wrote {len(rows)} rows -> {output_csv}")
    if missing_r2:
        print("WARNING: Missing R2 for the following R1 files:", file=sys.stderr)
        for name in missing_r2:
            print(f"  {name}", file=sys.stderr)

if __name__ == "__main__":
    main()
