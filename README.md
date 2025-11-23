# methyltna

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**methyltna** is a comprehensive bioinformatics pipeline for **TNA-EM-seq analysis** - a novel approach that combines **T**ranscriptome and methylome/**NA** analysis in a single sequencing experiment. The pipeline separates RNA and DNA fractions based on molecular barcodes, enabling simultaneous gene expression profiling and methylation analysis from the same sample.

### Key Features

- **🧬 TNA Deconvolution**: Separates RNA-barcoded and DNA-unbarcoded reads for dual analysis
- **📊 Gene Expression**: STAR alignment + FeatureCounts quantification for RNA reads
- **🔬 Methylation Analysis**: Biscuit EM-seq pipeline for DNA methylation calling
- **🧪 Variant Calling**: LoFreq low-frequency variant detection on RNA reads
- **📈 Comprehensive QC**: MultiQC integration with custom RNA barcode visualizations
- **💾 Smart Caching**: Reference genomes and indexes cached to save hours of computation
- **🔄 Paired-End Support**: Full paired-end read processing throughout pipeline

## Pipeline Summary

The pipeline performs the following analysis steps:

### 1. Input Processing
- **BCL Demultiplexing** ([`BCL-Convert`](https://support.illumina.com/sequencing/sequencing_software/bcl-convert.html) or [`bcl2fastq`](https://support.illumina.com/sequencing/sequencing_software/bcl2fastq-conversion-software.html)) - *Optional*
- **Samplesheet Validation** - nf-schema validation and channel creation

### 2. Read Quality Control & Trimming
- **Pre-trim QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- **Adapter Trimming** ([`TrimGalore`](https://github.com/FelixKrueger/TrimGalore))
  - Illumina adapter removal
  - Quality filtering (NextSeq 2-color chemistry support)
  - Minimum length filtering
- **Post-trim QC** ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))

### 3. RNA/DNA Separation (TNA Deconvolution)
- **RNA Barcode Extraction** ([`Cutadapt`](https://cutadapt.readthedocs.io/))
  - Extracts NNSR/mNNSR barcode sequences
  - **Barcoded reads** → RNA analysis path
  - **Unbarcoded reads** → DNA analysis path
- **Barcode Statistics** - MultiQC-compatible stacked bar visualizations

### 4. Reference Preparation (Smart Caching)
- **Download/Cache** - FASTA genome and GTF annotation files
- **STAR Index** ([`STAR`](https://github.com/alexdobin/STAR)) - Genome index for splice-aware alignment (~2+ hours, cached)
- **Biscuit Index** ([`Biscuit`](https://github.com/zwdzwd/biscuit)) - EM-seq index for methylation analysis
- **FASTA Index** ([`samtools faidx`](http://www.htslib.org/)) - Index for variant calling

### 5. RNA Analysis Path (Barcoded Reads)
- **Splice-Aware Alignment** ([`STAR`](https://github.com/alexdobin/STAR))
  - Uses GTF annotations for accurate splice junction detection
- **Format Conversion** ([`samtools view`](http://www.htslib.org/)) - SAM to BAM
- **Sorting** ([`samtools sort`](http://www.htslib.org/)) - Sort by coordinates
- **Gene Quantification** ([`FeatureCounts`](http://subread.sourceforge.net/))
  - Counts reads per gene
  - Produces counts matrix for differential expression
- **Variant Calling** ([`LoFreq`](https://csb5.github.io/lofreq/))
  - Low-frequency variant detection

### 6. DNA Methylation Analysis Path (Unbarcoded Reads)
- **EM-seq Alignment** ([`Biscuit`](https://github.com/zwdzwd/biscuit))
  - Optimized for enzymatic methylation detection
- **Duplicate Marking** ([`Picard MarkDuplicates`](https://broadinstitute.github.io/picard/))
  - Removes PCR duplicates
- **Alignment Statistics** ([`samtools`](http://www.htslib.org/))
  - Flagstat and idxstats reports
- **Methylation Calling** ([`Biscuit pileup`](https://github.com/zwdzwd/biscuit))
  - CpG, CHG, CHH methylation rates
  - VCF output format
- **Methylation QC** ([`Biscuit QC`](https://github.com/zwdzwd/biscuit))
  - Quality metrics (EM-seq specific)
- **Summary Statistics** - Aggregate methylation stats for MultiQC

### 7. Reporting
- **Comprehensive QC Report** ([`MultiQC`](http://multiqc.info/))
  - FastQC pre/post-trim metrics
  - TrimGalore adapter statistics
  - RNA barcode extraction rates (stacked bar plots)
  - STAR alignment statistics
  - FeatureCounts gene quantification
  - Biscuit methylation analysis
  - LoFreq variant statistics
  - Software versions

## Quick Start

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

### Test the Pipeline

```bash
nextflow run . -profile test,singularity --outdir results_test
```

## Usage

### Basic Example (Paired-End FASTQ)

```bash
nextflow run . \
    -profile singularity \
    --input samplesheet.csv \
    --genome_fasta /path/to/genome.fa \
    --annotation_gtf /path/to/annotation.gtf \
    --rna_barcode_config rna_barcodes.yaml \
    --outdir results
```

### Input Samplesheet

Create a CSV file with your paired-end FASTQ files:

**`samplesheet.csv`:**
```csv
sample,fastq_1,fastq_2
SAMPLE_1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
SAMPLE_2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

**Columns:**
- `sample`: Unique sample identifier
- `fastq_1`: Path to R1 reads (required)
- `fastq_2`: Path to R2 reads (required for paired-end)

### RNA Barcode Config

Create a YAML file defining your RNA barcode sequences:

**`rna_barcodes.yaml`:**
```yaml
adapter_name: "NNSR"
adapter_sequence: "NNSRAGTA"
reverse_complement_search: true
times: 1
```

### Core Parameters

#### Input/Output
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Path to FASTQ samplesheet (CSV) | `null` |
| `--bcl_input_dir` | Path to BCL directory (optional) | `null` |
| `--bcl_samplesheet` | Path to Illumina SampleSheet.csv | `null` |
| `--outdir` | Output directory | **Required** |

#### Reference Files
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--genome_fasta` | Path to genome FASTA (local or s3://) | **Required** |
| `--annotation_gtf` | Path to GTF annotation (local or s3://) | **Required** |
| `--reference_cache_dir` | Cache directory for references/indexes | `./references` |
| `--star_index` | Pre-built STAR index (optional) | `null` |
| `--biscuit_index` | Pre-built Biscuit index (optional) | `null` |

#### TNA Deconvolution
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rna_barcode_config` | YAML config with barcode sequences | **Required** |
| `--rna_min_overlap` | Minimum barcode overlap | `3` |
| `--rna_error_rate` | Barcode mismatch tolerance | `0.0` |
| `--rna_times` | Adapter search iterations | `1` |

#### Variant Calling (LoFreq)
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--lofreq_min_cov` | Minimum coverage depth | `10` |
| `--lofreq_min_bq` | Minimum base quality | `1` |
| `--lofreq_min_alt_bq` | Minimum alternate base quality | `1` |
| `--lofreq_sig` | Significance threshold | `0.01` |

### Advanced Options

#### BCL Demultiplexing

```bash
nextflow run . \
    -profile singularity \
    --bcl_input_dir /path/to/bcl \
    --bcl_samplesheet SampleSheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf annotation.gtf \
    --rna_barcode_config barcodes.yaml \
    --outdir results
```

#### Reference Caching Control

```bash
# Force re-download of references (ignores cache)
nextflow run . \
    -profile singularity \
    --input samplesheet.csv \
    --genome_fasta s3://bucket/genome.fa \
    --annotation_gtf s3://bucket/annotation.gtf \
    --force_redownload_references \
    --outdir results

# Force rebuild of indexes (ignores cached STAR/Biscuit indexes)
nextflow run . \
    -profile singularity \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf annotation.gtf \
    --force_rebuild_indexes \
    --outdir results
```

#### Use Pre-Built Indexes

```bash
nextflow run . \
    -profile singularity \
    --input samplesheet.csv \
    --genome_fasta genome.fa \
    --annotation_gtf annotation.gtf \
    --star_index /path/to/star/index \
    --biscuit_index /path/to/biscuit/index \
    --outdir results
```

### Execution Profiles

| Profile | Description |
|---------|-------------|
| `test` | Run with test data and resource limits |
| `docker` | Use Docker containers |
| `singularity` | Use Singularity containers (HPC) |
| `apptainer` | Use Apptainer containers |
| `conda` | Use Conda environments |

## Outputs

The pipeline generates comprehensive outputs organized by analysis type:

```
results/
├── fastqc/                          # FastQC reports
│   ├── raw/                         # Pre-trim QC
│   └── trimmed/                     # Post-trim QC
├── trimming/                        # TrimGalore outputs
│   └── *_trimming_report.txt
├── rna_deconvolution/               # RNA barcode extraction
│   ├── *_barcoded_R1.cutadapt.fastq # Barcoded R1
│   ├── *_barcoded_R2.cutadapt.fastq # Barcoded R2
│   ├── *_unbarcoded_R1.cutadapt.fastq
│   ├── *_unbarcoded_R2.cutadapt.fastq
│   ├── *.cutadapt.json              # JSON report
│   └── *.cutadapt.txt               # Text report
├── rna/                             # RNA barcode stats
│   └── rna_barcode_stats_mqc.txt
├── alignment/                       # Alignment outputs
│   ├── star/                        # RNA alignments
│   │   ├── *_barcoded.Aligned.out.sam
│   │   └── *_barcoded.Log.final.out
│   └── biscuit/                     # DNA alignments
│       └── *_unbarcoded.bam
├── gene_counts/                     # Gene quantification
│   ├── *_barcoded.featureCounts.tsv
│   └── *_barcoded.featureCounts.tsv.summary
├── samtools/                        # Alignment statistics
│   ├── *_barcoded.flagstat
│   ├── *_barcoded.idxstats
│   ├── *_unbarcoded.flagstat
│   └── *_unbarcoded.idxstats
├── methylation/                     # Methylation analysis
│   ├── *_unbarcoded.vcf.gz          # Methylation VCF
│   ├── *_unbarcoded_QC.txt          # Biscuit QC
│   └── biscuit_methylation_mqc.json # Summary stats
├── variants/                        # Variant calling
│   ├── *_barcoded.vcf.gz            # RNA variants
│   ├── *_barcoded.vcf.gz.tbi
│   └── lofreq_mqc.json              # Summary stats
├── multiqc/                         # Comprehensive QC
│   ├── multiqc_report.html          # Main report
│   └── multiqc_data/                # Underlying data
└── pipeline_info/                   # Execution info
    ├── execution_report.html
    ├── execution_timeline.html
    └── methyltna_software_mqc_versions.yml
```

### Key Output Files

#### Primary Outputs
- **`multiqc/multiqc_report.html`**: Comprehensive QC report
- **`gene_counts/*.featureCounts.tsv`**: Gene expression counts matrix
- **`methylation/*.vcf.gz`**: Methylation calls in VCF format
- **`variants/*.vcf.gz`**: Low-frequency variant calls

#### Intermediate Files
- **`rna_deconvolution/*_barcoded*.fastq`**: RNA-barcoded reads
- **`rna_deconvolution/*_unbarcoded*.fastq`**: DNA-unbarcoded reads
- **`alignment/star/*.sam`**: RNA alignments (splice-aware)
- **`alignment/biscuit/*.bam`**: DNA alignments (EM-seq)

#### QC Reports
- **`fastqc/raw/`**: Pre-trim FastQC reports
- **`fastqc/trimmed/`**: Post-trim FastQC reports
- **`trimming/*_trimming_report.txt`**: TrimGalore statistics
- **`rna/rna_barcode_stats_mqc.txt`**: RNA barcode extraction metrics
- **`samtools/*.flagstat`**: Alignment statistics
- **`methylation/*_QC.txt`**: Biscuit methylation QC

## Pipeline Architecture

For a detailed visualization of the pipeline workflow, see the [Pipeline Flowchart](docs/pipeline_flowchart.md).

**Key architectural features:**
- **Dual-path processing**: Separate RNA and DNA analysis pipelines
- **Smart caching**: Reference genomes and indexes cached across runs
- **Paired-end support**: Complete paired-end processing throughout
- **MultiQC integration**: Custom visualizations for RNA barcode statistics

### Workflow Overview

```
BCL/FASTQ → Read Trimming → TNA Deconvolution → ╔═══════════╗
                                                  ║ Barcoded  ║ → STAR → FeatureCounts → LoFreq
                                                  ╚═══════════╝
                                                  ╔═══════════╗
                                                  ║Unbarcoded ║ → Biscuit → Picard → Methylation
                                                  ╚═══════════╝
                                                         ↓
                                                     MultiQC
```

## Troubleshooting

### Common Issues

**Reference download failures (s3:// URLs):**
- Ensure AWS CLI authentication: `aws configure`
- Use `--force_redownload_references` to bypass cache
- Alternatively, download manually and provide local paths

**STAR index build memory errors:**
- STAR indexing requires ~32GB+ RAM for human genome
- Use pre-built indexes with `--star_index /path/to/index`
- Or use cloud/HPC resources with adequate memory

**RNA barcode extraction produces unexpected rates:**
- Verify barcode sequences in YAML config
- Check `--rna_error_rate` tolerance (default 0.0 is strict)
- Review cutadapt JSON reports for adapter matching

**Biscuit QC shows -nan values:**
- This is expected for EM-seq data (enzymatic vs. bisulfite conversion)
- Actual methylation data is in the BISCUIT_PILEUP VCF files
- Check VCF files for methylation call quality

### Resource Requirements

**Minimum recommended:**
- **CPUs**: 8+ cores
- **Memory**: 32GB+ (64GB+ for STAR indexing)
- **Storage**: 100GB+ for human genome analysis

**Typical usage (human genome, 10 samples):**
- **STAR alignment**: 32GB RAM, 8 CPUs, ~15 min/sample
- **Biscuit alignment**: 8GB RAM, 4 CPUs, ~30 min/sample
- **FeatureCounts**: 4GB RAM, 2 CPUs, ~5 min/sample
- **Total runtime**: ~2-4 hours (with parallelization)

## Credits

methyltna was originally written by [Marcus Viscardi](https://github.com/MarcusAtMotley) at [Motley Bio](https://motley.bio).

### Development Team
- **Marcus Viscardi** - Pipeline architecture, TNA deconvolution, MultiQC integration
- **Motley Bio** - Biological protocol development and validation

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For questions and support:
- [Open an issue](https://github.com/Motleybio-organization/modulestesting/issues) for bug reports
- Review [documentation](docs/pipeline_flowchart.md) for architecture details
- Check [troubleshooting](#troubleshooting) section for common problems

## Citations

If you use methyltna for your analysis, please cite the following tools:

### Core Tools
- **Nextflow**: Di Tommaso P, Chatzou M, Floden EW, et al. (2017). Nextflow enables reproducible computational workflows. Nat Biotechnol. 35(4):316-319. doi:10.1038/nbt.3820
- **MultiQC**: Ewels P, Magnusson M, Lundin S, Käller M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 32(19):3047-8. doi:10.1093/bioinformatics/btw354

### QC & Trimming
- **FastQC**: Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data. Available: http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- **Cutadapt**: Martin M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10-12. doi:10.14806/ej.17.1.200
- **TrimGalore**: Krueger F. TrimGalore: A wrapper around Cutadapt and FastQC. https://github.com/FelixKrueger/TrimGalore

### Alignment & Quantification
- **STAR**: Dobin A, Davis CA, Schlesinger F, et al. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 29(1):15-21. doi:10.1093/bioinformatics/bts635
- **Biscuit**: Zhou W, Dinh HQ, Ramjan Z, et al. (2018). DNA methylation loss in late-replicating domains is linked to mitotic cell division. Nat Genet. 50(4):591-602. doi:10.1038/s41588-018-0073-4
- **FeatureCounts**: Liao Y, Smyth GK, Shi W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 30(7):923-930. doi:10.1093/bioinformatics/btt656

### Variant Calling & Utilities
- **LoFreq**: Wilm A, Aw PP, Bertrand D, et al. (2012). LoFreq: a sequence-quality aware, ultra-sensitive variant caller. Nucleic Acids Research. 40(22):11189-11201. doi:10.1093/nar/gks918
- **Picard**: Broad Institute. Picard Toolkit. http://broadinstitute.github.io/picard/
- **SAMtools**: Li H, Handsaker B, Wysoker A, et al. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics. 25(16):2078-2079. doi:10.1093/bioinformatics/btp352

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

### Framework

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
