#!/usr/bin/env -S uv run --script
#./
# /// script
# requires-python = ">=3.12"
# dependencies = []
# ///
"""
Extract RNA barcode statistics from cutadapt JSON reports
for MultiQC custom content visualization.

This script parses JSON reports from RNA barcode extraction
and creates a TSV file that MultiQC can use to generate
a plot showing RNA barcode detection rates.
"""

import json
import argparse
import sys
from pathlib import Path


def parse_cutadapt_json(json_file):
    """Parse a cutadapt JSON file and extract RNA statistics."""
    try:
        with open(json_file, 'r') as f:
            data = json.load(f)

        # Extract basic info
        additional_info = data.get('additional_info', {})
        read_counts = data.get('read_counts', {})

        # Get sample name from input files or file name
        input_info = data.get('input', {})
        sample_name = Path(json_file).stem.replace('.cutadapt', '')

        # Extract key statistics
        total_reads = read_counts.get('input', 0)
        tagged_reads = additional_info.get('number_tagged', 0)
        double_tagged = additional_info.get('number_double_tagged', 0)
        single_tagged = additional_info.get('number_single_tagged', 0)
        untagged_reads = additional_info.get('number_untagged', 0)

        # Calculate percentages
        if total_reads > 0:
            pct_barcoded = (tagged_reads / total_reads) * 100
            pct_double_tagged = (double_tagged / total_reads) * 100
            pct_single_tagged = (single_tagged / total_reads) * 100
            pct_untagged = (untagged_reads / total_reads) * 100
        else:
            pct_barcoded = pct_double_tagged = pct_single_tagged = pct_untagged = 0

        return {
            'sample': sample_name,
            'total_reads': total_reads,
            'barcoded_reads': tagged_reads,
            'double_tagged': double_tagged,
            'single_tagged': single_tagged,
            'untagged_reads': untagged_reads,
            'pct_barcoded': pct_barcoded,
            'pct_double_tagged': pct_double_tagged,
            'pct_single_tagged': pct_single_tagged,
            'pct_untagged': pct_untagged
        }

    except Exception as e:
        print(f"Error parsing {json_file}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('json_files', nargs='+', help='Cutadapt JSON files to process')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')

    args = parser.parse_args()

    # Process all JSON files
    results = []
    for json_file in args.json_files:
        result = parse_cutadapt_json(json_file)
        if result:
            results.append(result)

    if not results:
        print("No valid JSON files processed", file=sys.stderr)
        sys.exit(1)

    # Write detailed TSV file
    with open(args.output, 'w') as f:
        f.write("Sample\tTotal Barcoded (%)\tDouble Tagged (%)\tSingle Tagged (%)\tUntagged (%)\n")
        for result in results:
            f.write(f"{result['sample']}\t{result['pct_barcoded']:.2f}\t{result['pct_double_tagged']:.2f}\t{result['pct_single_tagged']:.2f}\t{result['pct_untagged']:.2f}\n")

    # Create MultiQC custom content format
    mqc_output = args.output.replace('.tsv', '_mqc.txt')

    # Write in MultiQC custom content format as stacked bar plot with raw counts
    with open(mqc_output, 'w') as f:
        f.write("# plot_type: 'bargraph'\n")
        f.write("# section_name: 'RNA Barcode Detection'\n")
        f.write("# description: 'RNA barcode detection from TNA deconvolution showing barcoded vs unbarcoded reads'\n")
        f.write("# pconfig:\n")
        f.write("#   title: 'RNA Barcode Detection'\n")
        f.write("#   ylab: 'Number of reads'\n")
        f.write("#   cpswitch_counts_label: 'Number of reads'\n")
        f.write("#   cpswitch_c_active: False\n")
        f.write("#   stacking: 'normal'\n")
        f.write("Sample\tDouble_Tagged\tSingle_Tagged\tUnbarcoded\n")

        for result in results:
            f.write(f"{result['sample']}\t{result['double_tagged']}\t{result['single_tagged']}\t{result['untagged_reads']}\n")

    print(f"RNA statistics written to {args.output}")
    print(f"MultiQC data written to {mqc_output}")


if __name__ == '__main__':
    main()
