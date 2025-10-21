#!/usr/bin/env python3
"""
Extract methylation summary statistics from Biscuit VCF files for MultiQC.
Calculates average methylation levels for CpG, CHG, and CHH contexts.
"""

import sys
import os
import gzip
import json
from pathlib import Path
import argparse


def parse_biscuit_vcf(vcf_path):
    """Parse Biscuit VCF and extract methylation statistics."""
    stats = {
        'CpG': {'total': 0, 'methylated_sum': 0},
        'CHG': {'total': 0, 'methylated_sum': 0},
        'CHH': {'total': 0, 'methylated_sum': 0}
    }

    vcf_str = str(vcf_path)
    open_func = gzip.open if vcf_str.endswith('.gz') else open

    try:
        with open_func(vcf_str, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                # Extract context from INFO field
                info = fields[7]
                cx_match = [x for x in info.split(';') if x.startswith('CX=')]
                if not cx_match:
                    continue

                context = cx_match[0].split('=')[1]

                # Map Biscuit VCF context codes to display names
                # Biscuit uses: CG, CHG, CHH
                # Display uses: CpG, CHG, CHH (the 'p' represents phosphodiester bond)
                context_map = {'CG': 'CpG', 'CHG': 'CHG', 'CHH': 'CHH'}
                display_context = context_map.get(context, context)

                # Extract beta value (methylation level) from FORMAT field
                # BT field is the beta value (methylation level 0-1)
                format_fields = fields[8].split(':')
                sample_fields = fields[9].split(':')

                try:
                    bt_idx = format_fields.index('BT')
                    beta_value = float(sample_fields[bt_idx])

                    if display_context in stats:
                        stats[display_context]['total'] += 1
                        stats[display_context]['methylated_sum'] += beta_value
                except (ValueError, IndexError):
                    continue

    except Exception as e:
        print(f"Error parsing {vcf_path}: {e}", file=sys.stderr)
        return None

    # Calculate averages
    result = {}
    for context, data in stats.items():
        if data['total'] > 0:
            avg_methylation = (data['methylated_sum'] / data['total']) * 100
            result[context] = {
                'count': data['total'],
                'avg_methylation': round(avg_methylation, 2)
            }

    return result


def main():
    parser = argparse.ArgumentParser(description='Extract Biscuit methylation statistics for MultiQC')
    parser.add_argument('--input-dir', required=True, help='Directory containing Biscuit VCF files')
    parser.add_argument('--output', required=True, help='Output JSON file for MultiQC')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    vcf_files = list(input_dir.glob('*.vcf.gz'))

    if not vcf_files:
        print(f"Warning: No VCF files found in {input_dir}", file=sys.stderr)
        # Create empty output
        output_data = {
            'id': 'biscuit_methylation',
            'section_name': 'Biscuit Methylation Summary',
            'description': 'Average methylation levels across different cytosine contexts (EM-seq data)',
            'plot_type': 'table',
            'pconfig': {
                'id': 'biscuit_methylation_table',
                'title': 'Biscuit Methylation Summary'
            },
            'data': {}
        }
    else:
        # Process each VCF file
        all_stats = {}
        for vcf_file in vcf_files:
            sample_name = vcf_file.stem.replace('.vcf', '')
            stats = parse_biscuit_vcf(vcf_file)
            if stats:
                all_stats[sample_name] = stats

        # Create MultiQC custom content
        table_data = {}
        for sample, contexts in all_stats.items():
            table_data[sample] = {
                'CpG Sites': contexts.get('CpG', {}).get('count', 0),
                'CpG Methylation (%)': contexts.get('CpG', {}).get('avg_methylation', 0),
                'CHG Sites': contexts.get('CHG', {}).get('count', 0),
                'CHG Methylation (%)': contexts.get('CHG', {}).get('avg_methylation', 0),
                'CHH Sites': contexts.get('CHH', {}).get('count', 0),
                'CHH Methylation (%)': contexts.get('CHH', {}).get('avg_methylation', 0)
            }

        output_data = {
            'id': 'biscuit_methylation',
            'section_name': 'Biscuit Methylation Summary',
            'description': 'Average methylation levels across different cytosine contexts from Biscuit pileup VCF files',
            'plot_type': 'table',
            'pconfig': {
                'id': 'biscuit_methylation_table',
                'title': 'Biscuit Methylation Summary',
                'format': '{:,.0f}'
            },
            'data': table_data
        }

    # Write output
    with open(args.output, 'w') as f:
        json.dump(output_data, f, indent=2)

    print(f"Processed {len(vcf_files)} VCF files", file=sys.stderr)
    print(f"Generated MultiQC custom content: {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
