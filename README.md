# Motleybio-organization/modulestesting

[![GitHub Actions CI Status](https://github.com/Motleybio-organization/modulestesting/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Motleybio-organization/modulestesting/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Motleybio-organization/modulestesting/actions/workflows/linting.yml/badge.svg)](https://github.com/Motleybio-organization/modulestesting/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/Motleybio-organization/modulestesting)

## Introduction

**Motleybio-organization/modulestesting** is a comprehensive bioinformatics pipeline for **methylation sequencing analysis** with integrated RNA/DNA separation and gene quantification. Built on the nf-core framework, this pipeline implements a novel **hairpin-first architecture (v3.0)** that processes hairpin adapters before RNA/DNA separation, enabling cleaner downstream analysis and dual-aligner processing.

### Key Features

- **Hairpin-First Processing**: v2 protocol methylation chemistry with hairpin resolution before RNA/DNA separation
- **RNA/DNA Dual Processing**: Separates RNA-barcoded reads from DNA-unbarcoded reads for specialized analysis
- **Dual Aligner Strategy**: STAR (splice-aware, RNA) + Bowtie2 (genome, DNA)
- **Gene Quantification**: FeatureCounts for RNA expression analysis
- **Methylation Analysis**: Complete methylSNP workflow with bisulfite conversion
- **Variant Calling**: LoFreq for low-frequency variant detection
- **Comprehensive QC**: MultiQC reports integrating all analysis stages
- **Flexible Input**: Supports both FASTQ and raw BCL data

## Pipeline Summary

The pipeline performs the following analysis steps:

### 1. Input Processing
- **BCL Demultiplexing** ([`BCL-Convert`](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) or [`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)) - *Optional*
- **Samplesheet Validation** - nf-schema validation and channel creation
- **Pre-QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) - *Optional*

### 2. Hairpin Processing (v2 Protocol - Hairpin-First Architecture)
- **Illumina Adapter Trimming** ([`TrimGalore`](https://github.com/FelixKrueger/TrimGalore)) - Remove standard sequencing adapters
- **Hairpin Adapter Trimming** ([`TrimGalore`](https://github.com/FelixKrueger/TrimGalore)) - Remove custom hairpin adapters
- **Hairpin Resolution** (Python) - v2 protocol C→T conversion for methylation calling
- **Hairpin Statistics** - Generate MultiQC-compatible resolution metrics

### 3. RNA/DNA Separation (Optional)
- **RNA Barcode Extraction** ([`Cutadapt`](https://cutadapt.readthedocs.io/)) - Extract NNSR/mNNSR barcodes
  - Barcoded reads → RNA analysis path
  - Unbarcoded reads → DNA analysis path
- **RNA Statistics** - Generate barcode extraction metrics for MultiQC

### 4. Reference Preparation
- **Reference Download/Caching** - Smart caching for FASTA and GTF files
- **STAR Index** ([`STAR`](https://github.com/alexdobin/STAR)) - Genome index for splice-aware alignment
- **Bowtie2 Index** ([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) - Genome index for DNA alignment
- **FASTA Indexing** ([`samtools faidx`](http://www.htslib.org/)) - Index for variant calling

### 5. Alignment & Quantification

**Dual Processing Mode** (when RNA deconvolution enabled):
- **RNA Path (Barcoded Reads)**:
  - [`STAR`](https://github.com/alexdobin/STAR) - Splice-aware alignment with GTF annotations
  - [`samtools view`](http://www.htslib.org/) - SAM to BAM conversion
  - [`FeatureCounts`](http://subread.sourceforge.net/) - Gene-level quantification
- **DNA Path (Unbarcoded Reads)**:
  - [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - Genome alignment

**Single Processing Mode** (when RNA deconvolution skipped):
- [`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - All reads aligned to genome

### 6. Methylation Analysis (MethylSNP Processing)
- **Mark Unique Reads** (MarkUniread.py) - Identify uniquely mapping reads
- **Mark Duplicates** (MarkDup.py) - Remove PCR duplicates
- **Add Methylation Tags** (AddXMtag.py) - Add XM tags for methylation context
- **Bismark Extraction** ([`Bismark`](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)) - Extract CpG/CHG/CHH methylation rates

### 7. Variant Calling
- **BAM Indexing** ([`samtools index`](http://www.htslib.org/)) - Index processed BAM files
- **Variant Calling** ([`LoFreq`](https://csb5.github.io/lofreq/)) - Call low-frequency variants
- **Variant Statistics** - Aggregate VCF statistics for MultiQC

### 8. Reporting
- **Comprehensive QC Report** ([`MultiQC`](http://multiqc.info/)) - Integrated report with:
  - FastQC metrics
  - TrimGalore adapter statistics (Illumina + Hairpin)
  - Hairpin resolution metrics
  - RNA barcode extraction rates
  - STAR/Bowtie2 alignment statistics
  - FeatureCounts gene quantification summary
  - Bismark methylation analysis (CpG/CHG/CHH)
  - LoFreq variant calling statistics
  - Software versions

## Quick Start

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Test the Pipeline

```bash
nextflow run Motleybio-organization/modulestesting -profile test,docker --outdir results
```

## Usage

### Basic Example

```bash
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --input samplesheet.csv \\
    --genome_fasta /path/to/genome.fa \\
    --annotation_gtf /path/to/annotation.gtf \\
    --outdir results
```

### Input Options

#### Option 1: FASTQ Input (Recommended)

Prepare a samplesheet with your FASTQ files:

**`samplesheet.csv`:**
```csv
sample,fastq_1,fastq_2
SAMPLE_1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
SAMPLE_2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

**Columns:**
- `sample`: Unique sample identifier
- `fastq_1`: Path to forward reads (R1) - **Required**
- `fastq_2`: Path to reverse reads (R2) - **Required for paired-end**

#### Option 2: BCL Input (For raw Illumina data)

```bash
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --bcl_input_dir /path/to/bcl/directory \\
    --bcl_samplesheet /path/to/SampleSheet.csv \\
    --genome_fasta /path/to/genome.fa \\
    --annotation_gtf /path/to/annotation.gtf \\
    --outdir results
```

### Core Parameters

#### Input/Output
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to FASTQ samplesheet (CSV) | `null` |
| `--bcl_input_dir` | Path to BCL directory | `null` |
| `--bcl_samplesheet` | Path to Illumina SampleSheet.csv | `null` |
| `--outdir` | Output directory | **Required** |

#### Reference Files
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--genome_fasta` | Path to genome FASTA file (local or gs://) | **Required** |
| `--annotation_gtf` | Path to annotation GTF file (local or gs://) | **Required** |
| `--reference_cache_dir` | Directory for caching references and indexes | `./references` |
| `--star_index` | Pre-built STAR index directory | `null` |
| `--bowtie2_index` | Pre-built Bowtie2 index directory | `null` |

#### Processing Modes
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--skip_fastqc` | Skip FastQC on raw reads | `false` |
| `--skip_methylsnp_analysis` | Skip all methylSNP processing | `false` |
| `--skip_rna_deconvolution` | Skip RNA barcode extraction | `false` |
| `--rna_barcode_config` | Path to RNA barcode config YAML | `null` |

### Hairpin Processing Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--hairpin_illumina_adapter` | Illumina adapter sequence | Auto-detected |
| `--hairpin_adapter` | Custom hairpin adapter sequence | Project-specific |
| `--hairpin_min_overlap` | Minimum adapter overlap | `3` |
| `--hairpin_error_rate` | Adapter mismatch tolerance | `0.1` |

### RNA Barcode Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rna_barcode_config` | YAML config with barcode sequences | `null` |
| `--rna_min_overlap` | Minimum barcode overlap | `3` |
| `--rna_error_rate` | Barcode mismatch tolerance | `0.0` |
| `--rna_times` | Adapter search iterations | `1` |

**Example RNA barcode config (YAML):**
```yaml
adapter_name: "NNSR"
adapter_sequence: "NNSRAGTA"
reverse_complement_search: true
times: 1
```

### Variant Calling Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--lofreq_min_cov` | Minimum coverage depth | `10` |
| `--lofreq_min_bq` | Minimum base quality | `1` |
| `--lofreq_min_alt_bq` | Minimum alternate base quality | `1` |
| `--lofreq_sig` | Significance threshold | `0.01` |

### Advanced Options

#### Skip Specific Analyses

```bash
# Skip RNA deconvolution (process all reads as DNA)
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --input samplesheet.csv \\
    --skip_rna_deconvolution \\
    --genome_fasta genome.fa \\
    --annotation_gtf annotation.gtf \\
    --outdir results

# Skip methylSNP analysis entirely (FastQC + MultiQC only)
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --input samplesheet.csv \\
    --skip_methylsnp_analysis \\
    --outdir results
```

#### Reference Caching Control

```bash
# Force re-download of references (ignores cache)
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --input samplesheet.csv \\
    --genome_fasta gs://bucket/genome.fa \\
    --annotation_gtf gs://bucket/annotation.gtf \\
    --force_redownload_references \\
    --outdir results

# Force rebuild of indexes (ignores cached STAR/Bowtie2 indexes)
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    --input samplesheet.csv \\
    --genome_fasta genome.fa \\
    --annotation_gtf annotation.gtf \\
    --force_rebuild_indexes \\
    --outdir results
```

#### Custom Resource Limits

Create a `nf_limits.config` file:
```groovy
process {
  // Absolute resource ceilings
  resourceLimits = [ memory: 128.GB, cpus: 32, time: 168.h ]

  // Adjust specific processes
  withName: 'STAR_ALIGN' {
    memory = 64.GB
    cpus   = 16
  }

  withName: 'LOFREQ_CALLPARALLEL' {
    memory = 32.GB
    cpus   = 8
  }
}
```

Run with:
```bash
nextflow run Motleybio-organization/modulestesting \\
    -profile docker \\
    -c nf_limits.config \\
    --input samplesheet.csv \\
    --outdir results
```

### Execution Profiles

| Profile | Description |
|---------|-------------|
| `test` | Run with test data and resource limits |
| `docker` | Use Docker containers (recommended) |
| `singularity` | Use Singularity containers (for HPC) |
| `apptainer` | Use Apptainer containers |
| `conda` | Use Conda environments |

## Outputs

The pipeline generates comprehensive outputs organized by analysis type:

```
results/
├── fastqc/                          # FastQC reports (raw reads)
│   ├── sample1_fastqc.html
│   └── sample1_fastqc.zip
├── methylsnp_hairpin/               # Hairpin adapter trimming
│   ├── illumina_trimming/
│   │   └── sample1_trimming_report.txt
│   └── hairpin_trimming/
│       └── sample1_hairpin_trimming_report.txt
├── methylsnp/                       # Hairpin resolution + MethylSNP analysis
│   ├── final_sample1.Deconvolution.5mC              # Methylation report
│   ├── final_sample1_Deconvolution_R1.fq            # Hairpin-resolved reads
│   ├── sample1.hairpin_resolution_stats.txt         # Resolution statistics
│   ├── XM_sample1_barcoded.XMtag.sorted.bam         # Processed BAM with XM tags
│   ├── XM_sample1_barcoded.XMtag.sorted.bedGraph.gz # Methylation bedGraph
│   ├── XM_sample1_barcoded.XMtag.sorted.bismark.cov.gz # Bismark coverage
│   └── XM_sample1_barcoded.XMtag.sorted_splitting_report.txt # Bismark report
├── rna_deconvolution/               # RNA barcode extraction (if enabled)
│   └── cutadapt/
│       ├── sample1_barcoded.cutadapt.fastq      # Barcoded reads
│       ├── sample1_unbarcoded.cutadapt.fastq    # Unbarcoded reads
│       ├── sample1.cutadapt.json                # JSON report
│       └── sample1.cutadapt.txt                 # Text report
├── alignment/                       # Alignment outputs
│   ├── star/                        # RNA alignments (if dual mode)
│   │   ├── sample1_barcoded.Aligned.out.sam
│   │   └── sample1_barcoded.Log.final.out
│   └── bowtie2/                     # DNA alignments
│       ├── sample1_unbarcoded.sam
│       └── sample1_unbarcoded.bowtie2.log
├── gene_counts/                     # Gene quantification (if dual mode)
│   ├── sample1_barcoded.featureCounts.tsv
│   └── sample1_barcoded.featureCounts.tsv.summary
├── samtools/                        # SAMtools statistics
│   ├── sample1_barcoded.flagstat
│   ├── sample1_barcoded.idxstats
│   ├── sample1_unbarcoded.flagstat
│   └── sample1_unbarcoded.idxstats
├── variants/                        # Variant calling
│   ├── sample1_barcoded.vcf.gz
│   ├── sample1_barcoded.vcf.gz.tbi
│   └── sample1_unbarcoded.vcf.gz
├── hairpin/                         # Hairpin resolution MultiQC summaries
│   ├── hairpin_resolution_stats_mqc.yaml
│   └── hairpin_resolution_stats_mqc.tsv
├── rna/                             # RNA barcode MultiQC summaries
│   ├── rna_barcode_stats_mqc.txt
│   └── rna_barcode_stats.tsv
├── lofreq/                          # LoFreq MultiQC summaries
│   └── lofreq_mqc.json
├── multiqc/                         # Comprehensive QC report
│   ├── multiqc_report.html          # Main integrated report
│   └── multiqc_data/                # Underlying data
└── pipeline_info/                   # Pipeline execution info
    ├── execution_report.html        # Nextflow execution report
    ├── execution_timeline.html      # Execution timeline
    └── modulestesting_software_mqc_versions.yml  # Software versions
```

### Key Output Files

#### Primary Outputs
- **`multiqc/multiqc_report.html`**: Comprehensive QC report integrating all analysis stages
- **`methylsnp/*.bedGraph.gz`**: Methylation rates in bedGraph format
- **`methylsnp/*.bismark.cov.gz`**: Bismark coverage files (CpG, CHG, CHH methylation)
- **`variants/*.vcf.gz`**: Low-frequency variant calls
- **`gene_counts/*.featureCounts.tsv`**: Gene expression counts (RNA samples in dual mode)

#### Intermediate Files
- **`methylsnp/final_*_Deconvolution_R1.fq`**: Hairpin-resolved reads
- **`methylsnp/final_*.Deconvolution.5mC`**: Methylation reports from hairpin resolution
- **`rna_deconvolution/cutadapt/*_barcoded.cutadapt.fastq`**: RNA-barcoded reads
- **`rna_deconvolution/cutadapt/*_unbarcoded.cutadapt.fastq`**: DNA-unbarcoded reads
- **`alignment/star/*.sam`**: RNA alignments (splice-aware, dual mode only)
- **`alignment/bowtie2/*.sam`**: DNA alignments (genome)
- **`methylsnp/XM_*.XMtag.sorted.bam`**: Processed BAM files with XM methylation tags

#### Reports & Metrics
- **`methylsnp/*.hairpin_resolution_stats.txt`**: Per-sample hairpin resolution statistics
- **`hairpin/hairpin_resolution_stats_mqc.yaml`**: Aggregated hairpin stats for MultiQC
- **`rna_deconvolution/cutadapt/*.cutadapt.json`**: Barcode extraction metrics
- **`rna/rna_barcode_stats.tsv`**: Aggregated RNA barcode stats for MultiQC
- **`alignment/*/*.log`**: STAR and Bowtie2 alignment statistics
- **`samtools/*.flagstat`**: SAMtools flagstat reports for alignment QC
- **`methylsnp/*_splitting_report.txt`**: Bismark methylation extraction reports
- **`lofreq/lofreq_mqc.json`**: Aggregated variant calling statistics for MultiQC
- **`pipeline_info/`**: Detailed pipeline execution and version information

## Pipeline Architecture

For a detailed visualization of the pipeline workflow, see the [Pipeline Flowchart](docs/pipeline_flowchart.md).

**Key architectural features:**
- **Hairpin-First Processing**: v2 protocol processes hairpins before RNA/DNA separation
- **Conditional Branching**: Dual vs. single processing modes based on RNA deconvolution
- **Dual Aligner Strategy**: STAR (RNA, splice-aware) + Bowtie2 (DNA, genome)
- **Smart Caching**: Preserves expensive reference index builds across runs

## Configuration

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Example Parameter File

Create a `params.yml` file:
```yaml
# Input
input: 'samplesheet.csv'
outdir: 'results'

# References
genome_fasta: '/data/references/hg38.fa'
annotation_gtf: '/data/references/gencode.v44.annotation.gtf'
reference_cache_dir: './references'

# RNA deconvolution (optional)
rna_barcode_config: 'rna_barcodes.yaml'

# Processing options
skip_fastqc: false
skip_methylsnp_analysis: false

# Variant calling
lofreq_min_cov: 10
lofreq_sig: 0.01
```

Run with:
```bash
nextflow run Motleybio-organization/modulestesting -profile docker -params-file params.yml
```

## Troubleshooting

### Common Issues

**Reference download failures (gs:// URLs):**
- Ensure gcloud authentication is configured: `gcloud auth application-default login`
- Use `--force_redownload_references` to bypass cache
- Alternatively, download manually and provide local paths

**STAR index build memory errors:**
- STAR indexing requires ~32GB+ RAM for human genome
- Use pre-built indexes with `--star_index /path/to/star/index`
- Or use cloud/HPC resources with adequate memory

**MethylSNP processing failures:**
- Check that hairpin resolution produced valid output
- Ensure methylation reports (.5mC files) were generated
- Review individual sample logs in `work/` directory for specific errors

**LoFreq variant calling produces no variants:**
- Check BAM file has sufficient coverage (`--lofreq_min_cov`)
- Adjust significance threshold (`--lofreq_sig`) if too stringent
- Review samtools flagstat reports in MultiQC for alignment rates

**RNA deconvolution produces unexpected barcode rates:**
- Verify barcode sequences in config YAML
- Check `--rna_error_rate` tolerance (default 0.0 is strict)
- Review cutadapt JSON reports for adapter matching statistics

### Resource Requirements

**Minimum recommended resources:**
- **CPUs**: 8+ cores
- **Memory**: 32GB+ (64GB+ for STAR indexing)
- **Storage**: 100GB+ for human genome analysis

**Typical resource usage (human genome, 10 samples):**
- **STAR alignment**: 32GB RAM, 8 CPUs, ~15 min/sample
- **Bowtie2 alignment**: 8GB RAM, 4 CPUs, ~30 min/sample
- **LoFreq variant calling**: 16GB RAM, 4 CPUs, ~45 min/sample
- **Total runtime**: ~4-8 hours (depends on read depth and parallelization)

## Credits

Motleybio-organization/modulestesting was originally written by [Marcus Viscardi](https://github.com/MarcusAtMotley) at [Motley Bio](https://motley.bio).

### Development Team
- **Marcus Viscardi** - Pipeline architecture, methylSNP integration, hairpin-first design
- **Motley Bio** - Biological protocol development, validation, and testing

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For questions and support:
- [Open an issue](https://github.com/Motleybio-organization/modulestesting/issues) for bug reports or feature requests
- Review [documentation](docs/pipeline_flowchart.md) for pipeline architecture details
- Check [troubleshooting](#troubleshooting) section for common issues

## Citations

If you use Motleybio-organization/modulestesting for your analysis, please cite the following tools:

### Sequencing & QC
- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- **Cutadapt**: Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. doi:10.14806/ej.17.1.200
- **TrimGalore**: Krueger F. TrimGalore: A wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files. https://github.com/FelixKrueger/TrimGalore
- **MultiQC**: Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 32(19):3047-8. doi: 10.1093/bioinformatics/btw354

### Alignment
- **STAR**: Dobin A, Davis CA, Schlesinger F, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 29(1):15-21. doi:10.1093/bioinformatics/bts635
- **Bowtie2**: Langmead B, Salzberg SL. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods. 9(4):357-359. doi:10.1038/nmeth.1923

### Quantification & Analysis
- **FeatureCounts**: Liao Y, Smyth GK, Shi W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 30(7):923-930. doi:10.1093/bioinformatics/btt656
- **Bismark**: Krueger F, Andrews SR. (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics. 27(11):1571-1572. doi:10.1093/bioinformatics/btr167
- **LoFreq**: Wilm A, Aw PP, Bertrand D, et al. (2012). LoFreq: a sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Research. 40(22):11189-11201. doi:10.1093/nar/gks918

### Utilities
- **SAMtools**: Li H, Handsaker B, Wysoker A, et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics. 25(16):2078-2079. doi:10.1093/bioinformatics/btp352
- **BCL-Convert**: Illumina. BCL Convert v4.3.13. https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

### Framework

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
