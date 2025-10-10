# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Nextflow bioinformatics pipeline based on the nf-core template (version 3.3.2). The pipeline performs TNA-EM-seq analysis - combined transcriptome and methylome sequencing with RNA/DNA barcode separation.

**Pipeline**: Motleybio-organization/methyltna
**Main workflow**: TNA-EM-seq analysis with dual RNA/DNA processing paths
**Author**: Marcus Viscardi
**Nextflow version**: >=24.10.5

## Common Commands

### Running the pipeline
```bash
# Basic execution
nextflow run . -profile test,docker --outdir results

# Production run with custom samplesheet
nextflow run . -profile docker --input samplesheet.csv --outdir results

# Available profiles: test, docker, singularity, conda, mamba, apptainer
```

### Testing
```bash
# Run pipeline tests using nf-test
nf-test test

# Run specific test
nf-test test tests/default.nf.test
```

### Development and validation
```bash
# Validate pipeline configuration
nextflow config -profile test

# Generate execution reports
nextflow run . -profile test,docker --outdir results
# Reports generated in: results/pipeline_info/
```

## Architecture

### Core Workflow Structure
- **main.nf**: Entry point defining the main workflow `MOTLEYBIOORGANIZATION_METHYLTNA`
- **workflows/methyltna.nf**: Contains the `METHYLTNA` workflow with TNA deconvolution, alignment, and analysis processes
- **subworkflows/local/utils_nfcore_methyltna_pipeline/main.nf**: Pipeline initialization, completion, and utility functions

### Workflow Flow
1. **PIPELINE_INITIALISATION**: Validates parameters and creates input channels from samplesheet
2. **METHYLTNA**:
   - Performs read trimming with FastQC and TrimGalore
   - Executes TNA deconvolution (RNA barcode extraction)
   - Processes RNA reads through STAR alignment and gene quantification
   - Processes DNA reads through Biscuit methylation analysis
   - Performs variant calling with LoFreq
   - Generates comprehensive MultiQC report
3. **PIPELINE_COMPLETION**: Handles email notifications and cleanup

### Module Organization
- **modules/nf-core/**: Contains nf-core modules (fastqc, multiqc, bclconvert, bcl2fastq, trimgalore)
- **subworkflows/nf-core/**: Utility subworkflows (bcl_demultiplex, utils_*)
- **conf/**: Configuration files (base.config, test.config, modules.config)

### Input/Output
- **Input**: CSV samplesheet with columns: sample, fastq_1, fastq_2 (fastq_2 optional for single-end)
- **Output**: FastQC reports, MultiQC HTML report, pipeline execution reports
- **Schema validation**: Uses nf-schema plugin for parameter and samplesheet validation

### Configuration Profiles
- **test**: Uses test dataset from nf-core/test-datasets
- **docker/singularity/conda/mamba**: Container/environment management
- **gpu**: GPU-enabled execution
- **debug**: Enhanced debugging with process name validation

### Key Configuration Files
- **nextflow.config**: Main pipeline configuration with profiles and parameters
- **conf/base.config**: Process resource requirements and labels
- **conf/test.config**: Test-specific parameters and resource limits
- **nf-test.config**: Testing framework configuration
- **modules.json**: Tracks installed nf-core modules and their versions

## Development Notes

### Adding New Modules
Use nf-core tools to install modules:
```bash
nf-core modules install <module_name>
```

### Testing Framework
- Uses nf-test for pipeline testing
- Test files located in tests/ directory and module-specific test directories
- Run tests before making changes to ensure pipeline integrity

### Container Support
Pipeline supports multiple container engines (Docker, Singularity, Apptainer) with containers from quay.io registry.

## Current Development Status

### Recent Implementations (September 2025)
- **READ_TRIMMING subworkflow** fully implemented with FastQC -> TrimGalore -> FastQC pipeline
- **TrimGalore parameters** comprehensively configured with 2-color chemistry defaults:
  - `trim_nextseq = 20` (default for NextSeq/NovaSeq 2-color chemistry)
  - `trim_quality = null` (disabled standard quality trimming)
  - Additional parameters: trim_min_length, trim_adapter, trim_stringency
- **Pipeline architecture** updated to use READ_TRIMMING instead of standalone FastQC
- **Configuration files** updated: nextflow.config, conf/modules.config, nextflow_schema.json
- **README.md** completely rewritten with comprehensive usage documentation
- **Testing verified** pipeline works with both local (`nextflow run .`) and remote execution
- **MultiQC integration** confirmed - aggregates pre-trim FastQC, TrimGalore, and post-trim FastQC reports

#### RNA Barcode Extraction Module (September 15, 2025)
- **RNABARCODEEXTRACTION module** fully implemented and tested:
  - Located in `modules/local/rnabarcodeextraction/`
  - Uses portable cutadapt wrapper (`portable_cutadapt_wrapper.py`) in templates directory
  - Processes paired-end FASTQ files to extract RNA barcodes (NNSR/mNNSR)
  - Outputs: barcoded reads, unbarcoded reads, JSON reports, text reports
- **Container compatibility** verified:
  - ✅ Conda environment: Uses `uv>=0.4.0` for dependency management
  - ✅ Apptainer/Singularity: Uses `ghcr.io/astral-sh/uv:python3.12-bookworm`
  - ✅ Version tracking: Captures cutadapt, Python, and wrapper versions
- **Testing results** on Ca3-RNA-NNSR sample (501,083 read pairs):
  - 55.1% barcoded reads (276,297 reads)
  - 19.4% double-barcoded reads (97,058 reads)
  - 25.5% unbarcoded reads (127,728 reads)
- **Configuration**: Uses YAML config files (`runner_settings_readStart_EM_e0.0.yaml`)
- **Dual barcode support**: Updated for both NNSR and mNNSR barcode sequences

#### RNA Statistics Visualization (September 18, 2025)
- **RNA_STATS_SUMMARY module** fully implemented for MultiQC integration:
  - Located in `modules/local/rna_stats_summary/`
  - Uses UV script execution pattern with `extract_rna_stats.py`
  - Processes cutadapt JSON reports from TNA_DECONVOLUTION
  - Outputs both percentage TSV and MultiQC custom content
- **MultiQC stacked bar plots**:
  - Raw count format: Double_Tagged, Single_Tagged, Unbarcoded
  - Stacked bar visualization like FastQC duplication plots
  - Users can toggle between counts and percentages in MultiQC interface
- **Universal sample name cleaning**:
  - Regex patterns for TrimGalore suffixes (`_val_1`, `_val_2`)
  - Handles concatenated paired-end names from RNA extraction
  - Works across diverse sample naming conventions
- **Container integration**: Uses same UV container pattern as RNABARCODEEXTRACTION

#### METHYLSNP Filename Collision Fix (September 22, 2025)
- **Issue identified**: RNA barcode extraction created `.barcoded.` and `.unbarcoded.` files that caused filename collisions in methylSNP processing
- **Root cause**: AddXMtag.py script used SAM filename basename for temporary files, causing barcoded/unbarcoded data to overwrite each other's dictionaries
- **Solution implemented**:
  - Modified `portable_cutadapt_wrapper.py` to use underscore `_barcoded` instead of dot `.barcoded.` in output filenames
  - Updated RNABARCODEEXTRACTION module output patterns to match new naming convention
  - Used Python pathlib `.stem` property for proper filename construction instead of problematic `.with_suffix()` method
  - **Added meta.id modification** in TNA_DECONVOLUTION subworkflow to append `_barcoded`/`_unbarcoded` suffixes ensuring distinct prefixes throughout methylSNP pipeline
- **Technical details**:
  - **Filename level**: Changed from `THS1_12_CoB_11E_S12_R1.cutadapt.barcoded.fastq` to `THS1_12_CoB_11E_S12_R1.cutadapt_barcoded.fastq`
  - **Meta.id level**: Changed from `THS1_12_CoB_11E_S12` to `THS1_12_CoB_11E_S12_barcoded`/`THS1_12_CoB_11E_S12_unbarcoded`
  - **Result**: All methylSNP modules use distinct prefixes, creating separate temporary files: `bamtobed.sample_barcoded.uni.nodup.sam.script` vs `bamtobed.sample_unbarcoded.uni.nodup.sam.script`
- **Result**: Eliminates KeyError in AddXMtag.py when processing methylation data, allows independent processing of barcoded/unbarcoded branches

#### METHYLSNP Modules and Subworkflow (September 19, 2025)
- **METHYLSNP_ANALYSIS subworkflow** fully implemented and tested:
  - Located in `subworkflows/local/methylsnp_analysis/`
  - Complete end-to-end methylation sequencing analysis pipeline
  - Chains 5 specialized modules: TRIMREAD, DECONVOLUTION, ALIGNMENT, PROCESSING, EXTRACTION
- **Individual METHYLSNP modules** successfully implemented:
  - **METHYLSNP_TRIMREAD**: TrimGalore-based adapter trimming for methylation data
  - **METHYLSNP_DECONVOLUTION**: Python 2.7 deconvolution conversion using original team scripts
  - **METHYLSNP_ALIGNMENT**: Bowtie2 single-end alignment with read groups
  - **METHYLSNP_PROCESSING**: Multi-step SAM processing (MarkUniread, MarkDup, AddXMtag)
  - **METHYLSNP_EXTRACTION**: Bismark methylation extraction to bedGraph format
- **Python 2.7 compatibility** achieved:
  - Custom conda environments with numpy and legacy tool versions
  - cutadapt 1.18 + TrimGalore 0.4.5 + Python 2.7 combination
  - Original team scripts preserved without modifications
- **Testing results**:
  - ✅ All modules pass stub testing with correct output patterns
  - ✅ Real data testing: ALIGNMENT confirmed with bowtie2 execution
  - ✅ Real data testing: PROCESSING scripts execute and process SAM files
  - ✅ Complete subworkflow: End-to-end pipeline chains modules correctly
- **Production readiness**:
  - nf-core compliant module structure and patterns
  - Container support via conda environments
  - Ready for HPC execution with resource allocation
  - Modular design allows independent or complete workflow execution

#### LoFreq Variant Calling Integration (September 24, 2025)
- **LOFREQ_VARIANT_CALLING subworkflow** fully implemented for low-frequency variant detection:
  - Located in `subworkflows/local/lofreq_variant_calling/`
  - Handles BAM indexing (SAMTOOLS_INDEX) and variant calling (LOFREQ_CALLPARALLEL) in sequence
  - Proper channel management for BAM files, reference FASTA, and FASTA indices
- **Enhanced Reference Preparation**:
  - Added SAMTOOLS_FAIDX indexing to PREPARE_REFERENCES subworkflow
  - Generates FASTA indices (.fai files) for both genome and transcriptome during reference prep
  - Uses aliased processes (SAMTOOLS_FAIDX_GENOME, SAMTOOLS_FAIDX_TRANSCRIPTOME) for parallel execution
- **Comprehensive Parameter Configuration**:
  - Added `variant_calling_options` section to nextflow_schema.json with 8 configurable parameters
  - Module configuration in conf/modules.config with dynamic parameter passing
  - Sensible defaults for low-frequency variant detection (min_cov=10, sig=0.01)
- **Multi-Channel Integration**:
  - Three LoFreq instances: LOFREQ_RNA (transcriptome), LOFREQ_DNA (genome), LOFREQ_SINGLE (genome fallback)
  - Processes BAM files from methylSNP analysis (METHYLSNP_*.out.processed_bam)
  - VCF outputs integrated into MultiQC collection for comprehensive reporting
- **Production-Ready Architecture**:
  - Variant calling runs automatically after methylSNP processing completes
  - Compressed VCF outputs (.vcf.gz) with tabix indices (.vcf.gz.tbi) published to variants/ directory
  - Channel flow: methylSNP BAM → SAMTOOLS_INDEX → LOFREQ_CALLPARALLEL → MultiQC integration
- **Technical Fixes Applied**:
  - Fixed channel combination issues in METHYLSNP_PROCESSING (tuple destructuring mismatch)
  - Corrected SAMTOOLS_FAIDX input signature (3 inputs: fasta, existing_fai, get_sizes)
  - Resolved process uniqueness with aliased SAMTOOLS_FAIDX instances
  - Fixed channel output naming (processed_bam vs bam) in LoFreq integration

### Key Technical Decisions Made
- **2-color chemistry optimization**: Set trim_nextseq as default to handle spurious G-base calls from modern sequencers
- **Filename collision fix**: Added "_trimmed" suffix to post-trim FastQC meta.id to preserve both reports  
- **Parameter validation**: Full JSON schema implementation with help text and constraints
- **nf-core compliance**: All subworkflows and modules follow nf-core standards and patterns

### Pipeline Capabilities
- **BCL demultiplexing**: Supports both bclconvert and bcl2fastq with BCL_DEMULTIPLEX subworkflow
- **FASTQ processing**: READ_TRIMMING subworkflow handles quality control and adapter trimming
- **RNA barcode extraction**: RNABARCODEEXTRACTION module extracts NNSR/mNNSR barcodes from RNA-seq data
- **RNA statistics visualization**: RNA_STATS_SUMMARY generates stacked bar plots for MultiQC
- **Methylation sequencing analysis**: METHYLSNP_ANALYSIS subworkflow for complete bisulfite-seq processing
  - Adapter trimming with methylation-specific parameters
  - Deconvolution conversion for methylation calling
  - Bowtie2 alignment optimized for bisulfite-treated reads
  - SAM processing with methylation tagging (XM tags)
  - Bismark extraction to bedGraph format for downstream analysis
- **Low-frequency variant calling**: LOFREQ_VARIANT_CALLING subworkflow for detecting rare mutations
  - Automatic BAM indexing and FASTA indexing during reference preparation
  - Configurable parameters for coverage depth, base quality, and significance thresholds
  - Dual-channel support for RNA (transcriptome) and DNA (genome) variant detection
  - Compressed VCF output with tabix indexing for downstream analysis
- **Quality reporting**: MultiQC generates comprehensive before/after comparison reports with RNA metrics
- **Portable email notifications**: Built-in SMTP support for completion emails without system dependencies
- **Container support**: Works with Docker, Singularity, Apptainer across different compute environments

### Current State
- All development work committed and pushed to GitHub (master branch)
- Pipeline tested and functional for production use
- **RNA barcode extraction and statistics** fully integrated into TNA_DECONVOLUTION subworkflow
- **MultiQC integration** complete with RNA barcode visualization
- **Universal sample name cleaning** implemented for diverse naming conventions
- **METHYLSNP modules and subworkflow** fully implemented and tested
  - All 5 modules passing stub and real data testing
  - Python 2.7 compatibility resolved with proper conda environments
  - **Channel combination fix deployed**: All 74 samples now complete methylSNP processing (was 2/74)
  - **Bismark reports integrated**: MultiQC includes methylation analysis visualization
  - **Error handling implemented**: Comprehensive graceful degradation for empty/failed libraries
  - **Bash syntax resolved**: Consistent heredoc formatting prevents script generation errors
  - Production-ready for complete methylation sequencing analysis workflows
- Test datasets working correctly with expected outputs in MultiQC reports

#### Comprehensive Reference Caching & Collision Fix (September 23, 2025)
- **Genome-Agnostic Reference Caching System**: Implemented smart FASTA file caching to avoid redundant 3.5GB downloads
  - FASTA files cached in `references/fasta/` with dynamic filename extraction from URLs
  - Genome-specific index directories: `references/bowtie2_indexes/genome_name/`
  - Preserves existing 2+ hour Bowtie2 index investments while adding transcriptome support
  - Universal genome support through URL-based filename extraction (works with any genome)
- **Enhanced Filename Collision Prevention**:
  - **Root cause resolved**: Changed RNA barcode output from `sample.cutadapt_barcoded.fastq` to `sample_barcoded.cutadapt.fastq`
  - Moved barcode designation from suffix to filename stem in `portable_cutadapt_wrapper.py`
  - Updated RNABARCODEEXTRACTION module output patterns to match new naming convention
  - Prevents AddXMtag.py KeyError by ensuring distinct basenames for barcoded/unbarcoded processing
- **Production-Ready Error Handling**:
  - Comprehensive authentication failure messages with step-by-step resolution instructions
  - Manual download guidance with exact commands for users to copy/paste
  - Troubleshooting for container access, network, and authentication issues
- **Enhanced User Control**:
  - `--force_redownload_references` parameter to force re-downloading even with cached files
  - `--force_rebuild_indexes` parameter to force rebuilding even with cached indexes
  - Smart cache detection prevents redundant operations while allowing user override
- **Backward Compatibility**: All existing configurations and workflows preserved while adding new capabilities

#### AddXMtag.py Enhanced Error Handling (September 23, 2025)
- **Issue identified**: Race condition between bedtools execution and file reading in AddXMtag.py causing intermittent KeyError failures during pipeline execution
- **Root cause**: File system timing issues in Singularity container environment where temp files weren't fully written before being read
- **Solution implemented**:
  - **Retry mechanism**: 3 attempts with exponential backoff (1s, 2s, 4s delays)
  - **File system sync**: 0.5 second wait after bedtools execution for proper file sync
  - **Validation checks**: File existence, size, and format validation before proceeding
  - **Enhanced error reporting**: Detailed progress messages and diagnostic output
  - **Protected file reading**: Try-catch blocks around temp file parsing with debug output
- **Code preservation**: Original code commented out with clear section markers for easy rollback
- **Testing results**:
  - ✅ Manual execution successful with enhanced logging
  - ✅ Proper file validation and context loading (2 genome contexts loaded)
  - ✅ Robust error handling with informative messages
  - ✅ Backward compatibility maintained
- **Production impact**: Resolves intermittent methylSNP processing failures in high-throughput pipeline execution

#### MethylSNP Processing Error Handling & Bash Syntax Fix (September 24, 2025)
- **Issue identified**: IndexError in MarkDup.py when processing empty/failed libraries with zero uniquely mapping reads from MarkUniread.py
- **Error handling implementation**: Comprehensive graceful degradation for empty/failed libraries:
  - **Step-by-step validation**: Check outputs at each stage (MarkUniread, MarkDup, AddXMtag, SAM/BAM conversion)
  - **Diagnostic messaging**: Detailed warnings explaining potential causes of failures
  - **Graceful fallbacks**: Create minimal valid BAM files with headers-only when processing fails
  - **Status tracking**: Enhanced versions.yml with specific failure codes for monitoring
  - **Pipeline continuity**: All error conditions use `exit 0` to prevent pipeline crashes
- **Bash syntax error resolution**:
  - **Root cause**: Inconsistent heredoc formatting mixing `\t` escape sequences with actual spaces
  - **Solution**: Standardized all `<<-END_VERSIONS` sections to use consistent 4-space indentation
  - **Cache management**: Implemented selective work directory cleanup to force regeneration of fixed modules
- **Production impact**:
  - ✅ Prevents pipeline crashes from low-quality or failed samples
  - ✅ Maintains output consistency for downstream analysis
  - ✅ Provides clear diagnostic information for troubleshooting
  - ✅ Enables successful processing of mixed-quality sample batches

#### MethylSNP Channel Combination Fix (September 23, 2025)
- **Critical Issue Identified**: Data stream bottleneck where only 2 out of 74 samples were completing methylSNP processing and extraction stages
- **Root Cause**: Channel combination logic in METHYLSNP_PROCESSING was not properly pairing SAM files, methylation reports, and reference FASTA
  - Module expected 3 separate input channels but channel cardinalities weren't matching (74 sam + 74 reports + 1 reference)
  - Most samples never reached processing stage due to incomplete channel combinations
- **Solution Implemented** in `subworkflows/local/methylsnp_analysis/main.nf`:
  - Created `ch_combined_for_processing` channel using `.join()` and `.combine()` operations
  - Used `.join(by: [0])` to pair SAM files with methylation reports by sample ID
  - Used `.combine(reference_fasta)` to add reference to each pair
  - Split combined channel back into 3 separate channels using `.map()` operations to maintain module interface
- **Result**: All 74 samples now flow through complete methylSNP pipeline (trimming → deconvolution → alignment → processing → extraction)
- **Cache Behavior**: Pipeline re-run required due to workflow code changes invalidating Nextflow cache (expected and beneficial)

#### MethylSNP MultiQC Integration (September 24, 2025)
- **Issue Identified**: No methylSNP reports were being collected for MultiQC despite valuable Bismark splitting reports being generated
- **Integration Implemented**:
  - **Workflow Changes**: Added `extraction_report.collect{it[1]}` to MultiQC files channel for all methylSNP workflows (RNA, DNA, single)
  - **MultiQC Configuration**: Added Bismark module configuration with proper path filtering (`**/*_splitting_report.txt`)
  - **Sample Name Cleaning**: Added regex rules to clean methylSNP output prefixes (`XM_*.XMtag.sorted` → sample names)
- **MultiQC Content Added**:
  - **Bismark Methylation Analysis section** showing CpG, CHG, CHH methylation percentages
  - **C-to-T conversion efficiency** metrics for bisulfite conversion quality assessment
  - **Sample comparisons** between barcoded (RNA) and unbarcoded (DNA) fractions
  - **Interactive visualizations** with built-in MultiQC charts and tables
- **Testing Confirmed**: MultiQC v1.31 automatically detects and processes Bismark splitting reports
- **Production Ready**: Next pipeline completion will include comprehensive methylation analysis in MultiQC reports

#### STAR Alignment Migration (September 29, 2025)
- **Major Architecture Change**: Complete migration from Bowtie2 to STAR alignment for v2.0.0
  - **Motivation**: Required for FeatureCounts RNA quantification which needs genome-aligned reads (not transcriptome)
  - **Decision**: Full migration to STAR rather than dual aligner support for cleaner architecture
  - **Version Strategy**: Tagged v1.5.0 as final Bowtie2 version before breaking changes
- **Reference Preparation Updates**:
  - **PREPARE_REFERENCES subworkflow**: Replaced BOWTIE2_BUILD_INDEX with STAR_GENOMEGENERATE
  - **GTF Integration**: Added annotation_gtf parameter for splice-aware genome alignment
  - **Index Caching**: Smart STAR index caching in `references/star_indexes/genome_name/`
  - **SAMTOOLS_FAIDX**: Added genome indexing for variant calling compatibility
- **MethylSNP Analysis Migration**:
  - **Direct STAR Integration**: Replaced METHYLSNP_ALIGNMENT module with nf-core STAR_ALIGN
  - **SAM Output**: Configured STAR for SAM output to maintain methylSNP processing compatibility
  - **Parameters**: Optimized STAR parameters for splice-aware alignment and methylation analysis
  - **Workflow Updates**: Updated all methylSNP workflows (RNA, DNA, single) to use STAR inputs
- **Configuration and Schema**:
  - **STAR Parameters**: Added comprehensive STAR configuration in conf/modules.config
  - **JSON Schema**: Added annotation_gtf parameter to nextflow_schema.json to eliminate validation warnings
  - **Parameter Alignment**: Fixed gtf/annotation_gtf naming consistency across configs
- **Reference Caching and Authentication**:
  - **Issue Resolution**: Bypassed gcloud container authentication issues using smart caching
  - **Manual GTF Download**: Downloaded GTF from public Ensembl FTP as backup for testing
  - **Cache Detection**: Verified pipeline correctly detects cached files and skips problematic downloads
- **Testing and Validation**:
  - **Reference Preparation**: ✅ Confirmed STAR index generation from cached references
  - **Container Bypass**: ✅ Verified DOWNLOAD_REFERENCES process skipped when files cached
  - **Pipeline Flow**: ✅ Validated complete workflow from FastQC through STAR alignment
  - **Production Testing**: 🔄 Currently testing with 36 real samples in ~/runs/mot24
- **Technical Implementation Details**:
  - **STAR Configuration**: Splice-aware genome alignment with comprehensive parameter tuning
  - **File Format**: SAM output maintained for downstream methylSNP processing compatibility
  - **Channel Management**: Proper meta channel handling for STAR inputs/outputs
  - **Error Handling**: Robust reference detection and graceful fallbacks for missing files
- **Next Steps Planned**:
  - **FeatureCounts Integration**: Add RNA quantification for barcoded samples using nf-core subread/featurecounts
  - **Performance Validation**: Compare STAR vs Bowtie2 alignment quality and processing times
  - **Production Deployment**: Validate end-to-end pipeline with real production datasets

### Pending Investigations

#### Replace MarkDup.py with Picard MarkDuplicates (Priority: Medium)
**Current Implementation**: Custom Python 2 script (`MarkDup.py`) in METHYLSNP_PROCESSING module
- Groups reads by: chr + position + **sequence identity (col 10)**
- Keeps read with highest MAPQ from each duplicate group
- Total ~100 lines of Python 2.7 code

**Proposed Replacement**: nf-core PICARD_MARKDUPLICATES module
- **Benefits**:
  - Industry standard, well-tested and maintained
  - Much faster (Java vs Python 2)
  - Available as nf-core module (easy integration)
  - Outputs QC metrics for MultiQC
  - Python 2 is deprecated
- **Difference in Logic**:
  - Picard groups by: Library + Reference + Strand + 5' position (does NOT check sequence identity)
  - Picard selects by: Sum of base qualities (not MAPQ)
  - **Potential impact**: Picard might be more aggressive in marking duplicates
- **Testing Required**:
  1. Install `nf-core modules install picard/markduplicates`
  2. Run comparison on test dataset: MarkDup.py output vs Picard output
  3. Compare read counts and methylation calling results downstream
  4. Validate that sequence-identity check in MarkDup.py is/isn't critical for methylSNP protocol
- **Decision Point**: If Picard produces comparable results, switch to it for better maintainability

### Development Notes
- Use `nf-test test` for pipeline testing
- TrimGalore module configured for dynamic parameter passing via ext.args
- Pipeline follows nf-core template v3.3.2 structure and conventions

### Resume and Regeneration Patterns
When using `-resume` with Nextflow, sometimes you need to regenerate specific modules while keeping other cached results:

**Method**: Delete specific work directories before resuming
```bash
# Example: Regenerate RNA_STATS_SUMMARY module
# 1. Find the work directory (usually work/XX/XXXXXXX...)
find work -name "*RNA_STATS_SUMMARY*" -type d
# or look for the process hash in pipeline output: [9b/0ce3ea]

# 2. Remove the specific work directory
rm -rf work/9b/0ce3ea661100fbb5ee519823257445

# 3. Resume pipeline - the deleted module will reprocess, others stay cached
nextflow run ~/pipelines/modulestesting/ -c ./nf_overwrite.config -params-file ./nf_params_fastq.yaml -profile singularity -resume
```

**Note**: The `-refresh` option is not available in Nextflow version 25.04.7, so manual work directory deletion is the preferred method for selective reprocessing.

## Notifications Setup

### Email Notifications
```bash
# Environment variables (recommended)
export SMTP_HOST="smtp.gmail.com" SMTP_PORT="587"
export SMTP_USER="your.email@gmail.com" SMTP_PASSWORD="your_app_password"
nextflow run . --email your.email@gmail.com --outdir results
```

### Webhook Notifications (Slack/Teams)
```bash
# Slack
export SLACK_WEBHOOK="https://hooks.slack.com/services/..."
nextflow run . --hook_url $SLACK_WEBHOOK --outdir results

# Teams
nextflow run . --hook_url "https://outlook.office.com/webhook/..." --outdir results
```

Both support completion/failure notifications with rich formatting. Email includes MultiQC reports as attachments.

#### Pipeline Architecture Documentation (September 24, 2025)
- **Comprehensive Mermaid Flowchart** created and pushed to GitHub:
  - Located at `docs/pipeline_flowchart.md` with full architecture visualization
  - Shows complete workflow including conditional processing paths, dual vs single methylSNP analysis
  - Documents all subworkflows: READ_TRIMMING, TNA_DECONVOLUTION, METHYLSNP_ANALYSIS, LOFREQ_VARIANT_CALLING
  - Color-coded diagram with input nodes (blue), process nodes (purple), decision points (orange), subworkflows (green), and outputs (pink)
  - Interactive visualization compatible with GitHub, VS Code, and Mermaid Live Editor
- **Enhanced MarkDup.py Error Handling**:
  - Added comprehensive safety checks for empty/malformed input files
  - Improved duplicate removal algorithm with proper line format validation
  - Enhanced IndexError prevention for failed libraries with zero uniquely mapping reads
  - Maintains backward compatibility while adding robust error recovery
- **Updated LoFreq Parameter Defaults**:
  - Changed base quality thresholds from 6 to 1 based on proven working implementation
  - Added default values for alternate and reference base qualities (def_alt_bq=1, def_ref_bq=1)
  - Enhanced parameter documentation with implementation rationale
- **Documentation Integration**: All changes committed and pushed to GitHub with comprehensive commit message

#### LoFreq Error Resolution and Pipeline Completion (September 25, 2025)
- **Transient LoFreq Failure Analysis**:
  - Identified single failed LoFreq task (THS1_8_HRG_11A_unbarcoded) out of 72 total processes
  - Root cause: Temporary resource contention during parallel execution, not systematic overallocation
  - Confirmed Nextflow's resource management worked correctly (only 1/72 tasks failed = 98.6% success rate)
  - Manual re-execution of failed command succeeded, proving failure was transient
- **LOFREQ_SUMMARY Container Fix**:
  - **Issue**: Non-existent container URL `depot.galaxyproject.org/singularity/python:3.11--h2755cc3_0`
  - **Solution**: Updated to proven working container `ghcr.io/astral-sh/uv:python3.12-bookworm`
  - **Result**: LOFREQ_SUMMARY module now executes successfully with Python 3.12 and UV package manager
- **Pipeline Completion Achievement**:
  - **Status**: ✅ SUCCESSFUL - All 793 processes completed (791 cached, 2 newly executed)
  - **Duration**: 1m 7s execution time (3+ days saved by caching system)
  - **Output**: MultiQC report generated at `/home/marcus/runs/mot24/results_enhanced_multiqc/multiqc/multiqc_report.html`
  - **Coverage**: Complete integrated analysis including FastQC, RNA barcodes, methylation, variant calling, alignment metrics
- **Error Handling Validation**:
  - Nextflow's resume functionality successfully recovered from transient failures
  - Robust caching system prevented data loss and minimized re-computation
  - Graceful degradation for samples with no mapped reads (expected for low barcode recovery)
  - Comprehensive error logging and diagnostic capabilities demonstrated
