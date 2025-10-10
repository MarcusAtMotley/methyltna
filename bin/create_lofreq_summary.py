#!/usr/bin/env python3
"""
Generate LoFreq variant summary statistics for MultiQC.
This script parses VCF files from LoFreq and creates summary statistics.
"""

import sys
import os
import gzip
import json
from pathlib import Path
import argparse

def parse_vcf(vcf_file):
    """Parse VCF file and extract basic statistics."""
    stats = {
        'total_variants': 0,
        'snps': 0,
        'indels': 0,
        'transitions': 0,
        'transversions': 0,
        'mean_qual': 0.0,
        'mean_dp': 0.0
    }

    qual_sum = 0
    dp_sum = 0
    variant_count = 0

    try:
        if vcf_file.endswith('.gz'):
            file_handle = gzip.open(vcf_file, 'rt')
        else:
            file_handle = open(vcf_file, 'r')

        with file_handle as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 8:
                    continue

                ref = fields[3]
                alt = fields[4]
                qual = fields[5]
                info = fields[7]

                stats['total_variants'] += 1
                variant_count += 1

                # Classify variant type
                if len(ref) == 1 and len(alt) == 1:
                    stats['snps'] += 1
                    # Count transitions/transversions
                    if (ref == 'A' and alt == 'G') or (ref == 'G' and alt == 'A') or \
                       (ref == 'C' and alt == 'T') or (ref == 'T' and alt == 'C'):
                        stats['transitions'] += 1
                    else:
                        stats['transversions'] += 1
                else:
                    stats['indels'] += 1

                # Parse QUAL
                if qual != '.':
                    try:
                        qual_value = float(qual)
                        qual_sum += qual_value
                    except ValueError:
                        pass

                # Parse DP from INFO field
                for info_field in info.split(';'):
                    if info_field.startswith('DP='):
                        try:
                            dp_value = int(info_field.split('=')[1])
                            dp_sum += dp_value
                        except (ValueError, IndexError):
                            pass
                        break

        # Calculate means
        if variant_count > 0:
            stats['mean_qual'] = qual_sum / variant_count
            stats['mean_dp'] = dp_sum / variant_count

    except Exception as e:
        print(f"Error parsing {vcf_file}: {e}", file=sys.stderr)

    return stats

def create_multiqc_content(sample_stats):
    """Create MultiQC custom content from sample statistics."""

    # General stats table data
    general_stats = {}
    for sample, stats in sample_stats.items():
        general_stats[sample] = {
            'Total_Variants': stats['total_variants'],
            'SNPs': stats['snps'],
            'InDels': stats['indels'],
            'Ti_Tv_Ratio': round(stats['transitions'] / max(stats['transversions'], 1), 2),
            'Mean_QUAL': round(stats['mean_qual'], 1),
            'Mean_DP': round(stats['mean_dp'], 1)
        }

    # MultiQC custom content structure
    multiqc_data = {
        'id': 'lofreq_variants',
        'section_name': 'LoFreq Variant Calling',
        'description': 'Summary statistics from LoFreq low-frequency variant calling.',
        'plot_type': 'table',
        'pconfig': {
            'id': 'lofreq_stats_table',
            'title': 'LoFreq Variant Statistics',
            'save_file': True,
            'raw_data_fn': 'multiqc_lofreq_stats',
            'col1_header': 'Sample',
            'defaultsort': [{'idx': 1, 'dir': 'desc'}]  # Sort by Total_Variants descending
        },
        'data': general_stats,
        'headers': {
            'Total_Variants': {
                'title': 'Total Variants',
                'description': 'Total number of variants called',
                'scale': 'Greens',
                'format': '{:,.0f}'
            },
            'SNPs': {
                'title': 'SNPs',
                'description': 'Number of single nucleotide polymorphisms',
                'scale': 'Blues',
                'format': '{:,.0f}'
            },
            'InDels': {
                'title': 'InDels',
                'description': 'Number of insertions and deletions',
                'scale': 'Oranges',
                'format': '{:,.0f}'
            },
            'Ti_Tv_Ratio': {
                'title': 'Ti/Tv Ratio',
                'description': 'Transition/Transversion ratio',
                'scale': 'RdYlBu',
                'format': '{:.2f}'
            },
            'Mean_QUAL': {
                'title': 'Mean QUAL',
                'description': 'Mean variant quality score',
                'scale': 'viridis',
                'format': '{:.1f}'
            },
            'Mean_DP': {
                'title': 'Mean DP',
                'description': 'Mean depth of coverage',
                'scale': 'plasma',
                'format': '{:.1f}'
            }
        }
    }

    return multiqc_data

def main():
    parser = argparse.ArgumentParser(description='Generate LoFreq summary for MultiQC')
    parser.add_argument('--input-dir', '-i', required=True,
                        help='Directory containing LoFreq VCF files')
    parser.add_argument('--output', '-o', default='lofreq_mqc.json',
                        help='Output JSON file for MultiQC')

    args = parser.parse_args()

    input_path = Path(args.input_dir)
    sample_stats = {}

    # Find all VCF.gz files
    vcf_files = list(input_path.glob('**/*.vcf.gz'))

    if not vcf_files:
        print(f"No VCF.gz files found in {input_path}", file=sys.stderr)
        return 1

    print(f"Processing {len(vcf_files)} VCF files...")

    for vcf_file in vcf_files:
        # Extract sample name from filename
        sample_name = vcf_file.stem.replace('.vcf', '')
        print(f"Processing {sample_name}...")

        stats = parse_vcf(str(vcf_file))
        sample_stats[sample_name] = stats

    # Create MultiQC content
    multiqc_data = create_multiqc_content(sample_stats)

    # Write output
    with open(args.output, 'w') as f:
        json.dump(multiqc_data, f, indent=2)

    print(f"Created MultiQC data file: {args.output}")
    return 0

if __name__ == '__main__':
    sys.exit(main())