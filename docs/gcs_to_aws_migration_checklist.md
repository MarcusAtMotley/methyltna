# GCS → AWS Migration Checklist

**Date Created:** 2025-11-20
**Migration From:** GCS VM Development Environment
**Migration To:** AWS Infrastructure
**Pipeline:** methyltna v3.0+

---

## Executive Summary

This checklist guides the migration of the methyltna pipeline from the GCS VM development environment to AWS infrastructure. The VM contains ~6TB of data including pipeline code, reference files, and production run outputs. This document ensures all critical components are preserved or properly documented for regeneration on AWS.

---

## Pre-Migration Status (VM Audit Results)

### Git Repository Status
- ✅ **methyltna**: Clean working tree, all commits pushed to `origin/master`
- ✅ **stemloopTNA-pipeline** (formerly modulestesting): All changes committed and pushed (commit fe387eb)

### Data Inventory
- **Pipeline caches**: 3.7GB (Nextflow work directory - regenerable)
- **Reference files**: 44GB (genome, indexes, annotations)
- **Production runs**: ~5.7TB (mot26: 4.4TB, mot30: 821GB, others: ~500GB)
- **Credentials**: Found but not in git (will need to recreate on AWS)

---

## Phase 1: Pre-Closure Actions (ON VM)

### 1.1 Handle stemloopTNA-pipeline Repository ✅ COMPLETED

**Location:** `/home/marcus/pipelines/stemloopTNA-pipeline` (formerly `modulestesting`)
**GitHub:** `https://github.com/Motleybio-organization/stemloopTNA-pipeline`

- [x] Review uncommitted changes ✅
- [x] **Committed and pushed all changes** (commit fe387eb) ✅
  - Cloud reference caching system (DOWNLOAD_CLOUD_CACHE module)
  - Custom container management with build scripts
  - Enhanced PREPARE_REFERENCES subworkflow
  - samtools/view nf-core module
  - Documentation updates (CLAUDE.md, README.md)
  - 23 files changed, 2412 insertions

- [x] **Repository renamed** from `modulestesting` to `stemloopTNA-pipeline` ✅
  - Local directory updated: `~/pipelines/stemloopTNA-pipeline`
  - Remote URL updated: `https://github.com/Motleybio-organization/stemloopTNA-pipeline.git`

- [x] **Container files properly managed** ✅
  - Container recipe (`.def`) and build scripts tracked in git
  - Binary (`.sif`) properly gitignored
  - Documentation added to `containers/README.md`

### 1.2 Preserve Run Outputs (CRITICAL DECISION)

**Location:** `~/runs/` (5.7TB total)

**Question:** Which run outputs are irreplaceable?

- [ ] **Option A: Archive all runs (expensive)**
  - Cost: ~$150/month S3 Standard or ~$30/month S3 Glacier
  - Command:
    ```bash
    # Upload to GCS first (while authenticated)
    gsutil -m cp -r ~/runs/ gs://motleybio-archive/vm-runs/
    # Then transfer to S3 later
    aws s3 sync gs://motleybio-archive/vm-runs/ s3://motleybio-archive/vm-runs/
    ```

- [ ] **Option B: Archive selected outputs only (recommended)**
  - [ ] Keep MultiQC reports (small, valuable):
    ```bash
    find ~/runs -name "multiqc_report.html" -exec cp --parents {} /tmp/run_reports/ \;
    gsutil -m cp -r /tmp/run_reports/ gs://motleybio-archive/
    ```
  - [ ] Keep final BAMs from critical runs (list below):
    - [ ] mot26 BAMs? (Y/N): _____
    - [ ] mot30 BAMs? (Y/N): _____
    - [ ] mot27 BAMs? (Y/N): _____

  - [ ] Document run parameters for regeneration:
    ```bash
    find ~/runs -name ".nextflow.log" -o -name "params.json" -exec cp --parents {} /tmp/run_params/ \;
    ```

- [ ] **Option C: Delete all (regenerate on AWS)**
  - [ ] Confirmed: All runs can be regenerated from source data? (Y/N): _____
  - [ ] Backup only samplesheets:
    ```bash
    find ~/runs -name "samplesheet*.csv" -exec cp --parents {} /tmp/samplesheets/ \;
    ```

### 1.3 Transfer Reference Files to S3

**Current sources (GCS):**
- `gs://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa` (3.1GB)
- `gs://motleybio/Resources/RSEM_Genome/RSEM_hg38.transcripts.fa`
- `gs://motleybio/Resources/GTF/Homo_sapiens.GRCh38.112.chr_label.gtf` (1.5GB)

**Options:**

- [ ] **Option A: Transfer pre-built indexes (44GB, saves time)**
  ```bash
  # Create S3 bucket structure
  aws s3 mb s3://motleybio-references

  # Transfer cached references from VM
  aws s3 sync ~/pipelines/methyltna/references/ s3://motleybio-references/ \
    --exclude "*.log" --exclude ".nextflow*"

  # Estimated time: 30-60 minutes
  # Estimated cost: ~$5 transfer + ~$1/month storage
  ```

- [ ] **Option B: Transfer only FASTA/GTF (4.6GB, rebuild indexes on AWS)**
  ```bash
  # Transfer source files only
  aws s3 sync ~/pipelines/methyltna/references/fasta/ s3://motleybio-references/fasta/
  aws s3 cp ~/pipelines/methyltna/references/Homo_sapiens.GRCh38.112.chr_label.gtf \
    s3://motleybio-references/gtf/

  # Let pipeline rebuild indexes on first AWS run (~3 hours, ~$2 compute)
  ```

- [ ] **Option C: Use public sources (no transfer, download from Ensembl)**
  - Pipeline will download from Ensembl on first run
  - May be slower depending on network

**Selected Option:** _____ (A/B/C)

### 1.4 Backup Custom Scripts and Analysis

- [ ] Search for custom scripts not in git:
  ```bash
  find ~ -type f \( -name "*.py" -o -name "*.R" -o -name "*.sh" -o -name "*.ipynb" \) \
    -not -path "*/.*" -not -path "*/work/*" -not -path "*/miniconda3/*" \
    -not -path "*/miniforge3/*" -mtime -90 > /tmp/custom_scripts_list.txt

  cat /tmp/custom_scripts_list.txt
  ```

- [ ] Review list and backup any important analysis scripts:
  ```bash
  # Download to local machine
  scp marcus@vm-name:/path/to/script.py ~/local/backup/
  ```

### 1.5 Document Current Environment

- [ ] Export Nextflow configuration:
  ```bash
  cd ~/pipelines/methyltna
  nextflow config -profile singularity > /tmp/current_vm_config.txt
  ```

- [ ] Document VM specifications (for AWS sizing):
  ```bash
  lscpu | grep -E "^CPU\(s\)|Model name|Thread|Core"
  free -h
  df -h
  # Record: CPUs: ____ RAM: ____ Storage: ____
  ```

- [ ] List installed software versions:
  ```bash
  nextflow -version > /tmp/software_versions.txt
  singularity --version >> /tmp/software_versions.txt
  docker --version >> /tmp/software_versions.txt 2>&1
  conda --version >> /tmp/software_versions.txt 2>&1
  ```

### 1.6 Final Git Verification

- [ ] Verify methyltna repository:
  ```bash
  cd ~/pipelines/methyltna
  git status
  git log --oneline -5
  git log --branches --not --remotes --oneline  # Should be empty
  ```

- [x] Verify stemloopTNA-pipeline (after handling changes in 1.1): ✅
  ```bash
  cd ~/pipelines/stemloopTNA-pipeline
  git status  # Clean - only .claude/settings.local.json (local config, not tracked)
  git log --branches --not --remotes --oneline  # Empty - all pushed
  ```

- [ ] Clone to local machine as final backup:
  ```bash
  # On local machine
  git clone https://github.com/MarcusAtMotley/methyltna.git
  git clone https://github.com/Motleybio-organization/stemloopTNA-pipeline.git
  ```

---

## Phase 2: AWS Setup (AFTER VM Closure)

### 2.1 AWS Account Configuration

- [ ] Create/verify AWS account access
- [ ] Set up AWS credentials locally:
  ```bash
  aws configure
  # AWS Access Key ID: __________
  # AWS Secret Access Key: __________
  # Default region: us-east-1 (or preferred)
  # Default output format: json
  ```

- [ ] Test AWS access:
  ```bash
  aws sts get-caller-identity
  aws s3 ls
  ```

### 2.2 Create S3 Infrastructure

- [ ] Create S3 buckets:
  ```bash
  # Reference files
  aws s3 mb s3://motleybio-references --region us-east-1

  # Pipeline results
  aws s3 mb s3://motleybio-results --region us-east-1

  # Archive (if needed)
  aws s3 mb s3://motleybio-archive --region us-east-1
  ```

- [ ] Set up bucket policies and encryption:
  ```bash
  # Enable versioning on references bucket
  aws s3api put-bucket-versioning \
    --bucket motleybio-references \
    --versioning-configuration Status=Enabled

  # Enable encryption
  aws s3api put-bucket-encryption \
    --bucket motleybio-references \
    --server-side-encryption-configuration \
    '{"Rules":[{"ApplyServerSideEncryptionByDefault":{"SSEAlgorithm":"AES256"}}]}'
  ```

### 2.3 AWS Batch Setup

**Documentation:** See `docs/aws_batch_prerequisites.md` and `docs/aws_batch_setup.md`

- [ ] Create VPC and subnets (or use existing)
- [ ] Create AWS Batch compute environment
- [ ] Create job queue
- [ ] Create IAM roles:
  - [ ] Batch service role
  - [ ] ECS task execution role
  - [ ] EC2 instance role (with S3 access)

- [ ] Test AWS Batch:
  ```bash
  aws batch describe-compute-environments
  aws batch describe-job-queues
  ```

### 2.4 Update Pipeline for AWS

- [ ] Clone repository on development machine:
  ```bash
  git clone https://github.com/MarcusAtMotley/methyltna.git
  cd methyltna
  git checkout -b aws-migration
  ```

- [ ] Update `nextflow.config` with S3 paths:
  ```groovy
  params {
      // OLD (GCS):
      // genome_fasta = "gs://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

      // NEW (S3):
      genome_fasta = "s3://motleybio-references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa"
      transcriptome_fasta = "s3://motleybio-references/fasta/RSEM_hg38.transcripts.fa"
      annotation_gtf = "s3://motleybio-references/gtf/Homo_sapiens.GRCh38.112.chr_label.gtf"

      // Update cache directory to use S3
      reference_cache_dir = "s3://motleybio-references"
  }
  ```

- [ ] Update `modules/local/download_references/main.nf`:
  - [ ] Change container from `google/cloud-sdk:alpine` to `amazon/aws-cli:latest`
  - [ ] Replace `gsutil cp` with `aws s3 cp`
  - [ ] Update error messages for AWS context

- [ ] Update AWS Batch profile in `conf/awsbatch.config`:
  - [ ] Set work directory: `workDir = 's3://motleybio-results/work'`
  - [ ] Configure job queue name
  - [ ] Set resource limits

- [ ] Commit and push changes:
  ```bash
  git add -A
  git commit -m "feat: migrate from GCS to AWS S3 storage"
  git push origin aws-migration
  ```

- [ ] Create pull request and merge after review

### 2.5 First AWS Test Run

- [ ] Run test profile to verify setup:
  ```bash
  nextflow run MarcusAtMotley/methyltna \
    -profile test,awsbatch \
    --outdir s3://motleybio-results/test-run-001/ \
    -work-dir s3://motleybio-results/work/ \
    -resume
  ```

- [ ] Monitor execution:
  ```bash
  # Watch AWS Batch console
  # Check CloudWatch logs
  # Monitor S3 for outputs
  ```

- [ ] Verify outputs:
  ```bash
  aws s3 ls s3://motleybio-results/test-run-001/multiqc/
  aws s3 cp s3://motleybio-results/test-run-001/multiqc/multiqc_report.html .
  ```

- [ ] Document first-run performance:
  - Total runtime: _____ minutes
  - Cost estimate: $_____
  - Index generation time: _____ minutes
  - Issues encountered: _____________________

### 2.6 Production Run Validation

- [ ] Prepare production samplesheet:
  ```bash
  # Upload to S3
  aws s3 cp samplesheet.csv s3://motleybio-results/samplesheets/
  ```

- [ ] Launch production run:
  ```bash
  nextflow run MarcusAtMotley/methyltna \
    -r master \
    -profile production,awsbatch \
    --input s3://motleybio-results/samplesheets/samplesheet.csv \
    --outdir s3://motleybio-results/production-run-001/ \
    -work-dir s3://motleybio-results/work/ \
    -resume
  ```

- [ ] Validate against VM results:
  - [ ] Compare MultiQC metrics
  - [ ] Spot-check BAM file sizes/checksums
  - [ ] Verify variant calling results
  - [ ] Check methylation analysis outputs

---

## Phase 3: Decommissioning VM

### 3.1 Final Verification

- [ ] Confirm all critical data backed up:
  - [x] methyltna repository pushed to GitHub ✅
  - [x] stemloopTNA-pipeline repository pushed to GitHub ✅
  - [ ] Reference files in S3
  - [ ] Critical run outputs archived (if applicable)
  - [ ] Custom scripts downloaded

- [ ] Confirm AWS pipeline working:
  - [ ] Test run completed successfully
  - [ ] Production run validated
  - [ ] Reference indexes cached in S3

### 3.2 Clean Shutdown

- [ ] Stop any running Nextflow processes:
  ```bash
  ps aux | grep nextflow
  # Kill if needed: kill <PID>
  ```

- [ ] Export final logs:
  ```bash
  tar -czf vm-logs-$(date +%Y%m%d).tar.gz ~/.nextflow.log ~/runs/*/.nextflow.log
  # Upload to S3
  aws s3 cp vm-logs-*.tar.gz s3://motleybio-archive/vm-final-backup/
  ```

- [ ] Document VM configuration for reference:
  ```bash
  # Save to file and upload
  echo "VM Hostname: $(hostname)" > vm-config.txt
  echo "VM IP: $(hostname -I)" >> vm-config.txt
  echo "OS: $(cat /etc/os-release | grep PRETTY_NAME)" >> vm-config.txt
  echo "Nextflow: $(nextflow -version 2>&1 | head -1)" >> vm-config.txt
  aws s3 cp vm-config.txt s3://motleybio-archive/vm-final-backup/
  ```

### 3.3 VM Deletion

- [ ] Final confirmation checklist:
  - [ ] All code in GitHub? ✓
  - [ ] All reference files accessible on AWS? ✓
  - [ ] AWS pipeline tested and validated? ✓
  - [ ] No irreplaceable data left on VM? ✓

- [ ] Delete VM instance:
  ```bash
  # GCP command (adjust for your infrastructure)
  gcloud compute instances delete <vm-name> --zone=<zone>
  ```

---

## Phase 4: Post-Migration (Optional)

### 4.1 Cost Optimization

- [ ] Monitor first month AWS costs
- [ ] Set up CloudWatch billing alarms
- [ ] Review and optimize:
  - [ ] Spot instances for compute environments
  - [ ] S3 lifecycle policies (move old runs to Glacier)
  - [ ] Reserved capacity if usage is predictable

### 4.2 Documentation Updates

- [ ] Update README.md with AWS-specific instructions
- [ ] Update CLAUDE.md with AWS context
- [ ] Create `docs/aws_best_practices.md` with learnings

### 4.3 Team Communication

- [ ] Notify team of migration completion
- [ ] Share new AWS execution commands
- [ ] Update any CI/CD pipelines
- [ ] Archive VM-specific documentation

---

## Reference Information

### Current GCS Paths (for reference)
```
genome_fasta:        gs://motleybio/Resources/reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
transcriptome_fasta: gs://motleybio/Resources/RSEM_Genome/RSEM_hg38.transcripts.fa
annotation_gtf:      gs://motleybio/Resources/GTF/Homo_sapiens.GRCh38.112.chr_label.gtf
```

### New AWS Paths (to be created)
```
genome_fasta:        s3://motleybio-references/fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa
transcriptome_fasta: s3://motleybio-references/fasta/RSEM_hg38.transcripts.fa
annotation_gtf:      s3://motleybio-references/gtf/Homo_sapiens.GRCh38.112.chr_label.gtf
```

### Cost Estimates
- **Reference storage (44GB)**: ~$1/month (S3 Standard)
- **First pipeline run** (with index generation): ~$2-5 compute
- **Subsequent runs**: Variable based on sample count and size
- **Data transfer (GCS→S3)**: ~$0.12/GB egress from GCS

### Key Contacts
- Pipeline Developer: Marcus Viscardi
- GitHub Repository: https://github.com/MarcusAtMotley/methyltna
- Stem-Loop TNA Pipeline: https://github.com/Motleybio-organization/stemloopTNA-pipeline

---

## Completion Sign-Off

- [ ] Phase 1 (Pre-Closure) completed by: __________ Date: __________
- [ ] Phase 2 (AWS Setup) completed by: __________ Date: __________
- [ ] Phase 3 (Decommissioning) completed by: __________ Date: __________
- [ ] VM successfully deleted: __________ Date: __________

**Migration Status:** ☐ In Progress ☐ Complete

**Notes:**
```
[Add any migration-specific notes, issues encountered, or deviations from plan here]
```
