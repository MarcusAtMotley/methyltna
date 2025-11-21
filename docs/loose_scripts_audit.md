# Loose Scripts and Projects Audit

**Date:** 2025-11-21
**Purpose:** Identify any important code/data not tracked in git repositories before VM closure

---

## ✅ Already Backed Up in GitHub

### 1. ~/pipelines/methyltna
- **Repository:** https://github.com/MarcusAtMotley/methyltna
- **Status:** ✅ Clean, all commits pushed
- **Latest commit:** 4f611d7 - "fix: correct S3 bucket name in migration checklist"

### 2. ~/pipelines/stemloopTNA-pipeline
- **Repository:** https://github.com/Motleybio-organization/stemloopTNA-pipeline
- **Status:** ✅ Clean, all commits pushed
- **Latest commit:** fe387eb - "feat: add cloud reference caching and container management"

### 3. ~/projects/BarcodeExtracting
- **Repository:** https://github.com/MarcusAtMotley/BarcodeExtracting
- **Status:** ✅ Clean, all commits pushed
- **Latest commit:** 358dfb5 - "changed the main NNSR_EMaware_setting.yaml to better settings"
- **Contains:**
  - paired_fastq_finder.py / paired_fastq_finder_gcloud.py
  - portable_cutadapt_wrapper.py
  - portable_cutadapt_multirunner.py
  - Jupyter notebooks for testing and visualization

### 4. ~/projects/hairpinQuickAligner
- **Repository:** https://github.com/MarcusAtMotley/hairpinQuickAligner
- **Status:** ✅ Clean, all commits pushed
- **Latest commit:** 151473d - "init"
- **Contains:**
  - portable_hairpinAligner.py
  - HairpinAlignment.ipynb and related notebooks
  - FASTQ processing utilities

### 5. ~/dotfiles
- **Repository:** https://github.com/MarcusAtMotley/dotfiles
- **Status:** ✅ Clean, all commits pushed
- **Contains:**
  - .bashrc
  - .bash_aliases
  - Personal shell configuration

---

## ❌ NOT Backed Up - Requires Action

### ⚠️ ~/sunil_offload/ - Machine Learning Deconvolution Project

**Status:** NOT a git repository
**Last Modified:** 2025-11-12 (recent!)
**Total Size:** ~29MB code + ~156MB outputs
**Priority:** HIGH - Contains trained models and analysis code

#### Contents:

**Python Scripts:**
- `Motley_FM_Model_and_Cell_Type_Deconvolution.py` (42KB, 1,375 lines)
  - Cross-modal RNA-methylation foundation model
  - TCGA data processing with RNA expression and methylation
  - Transformer-based architecture (d_model=512, 8 heads, 6 layers)
  - Numeric ID prediction from sample names
  - Test-time augmentation (TTA) with 4-view averaging

**Training Configuration:**
- `Motley_Train_Model.sh` (1KB)
  - Batch size: 32
  - Learning rate: 2e-4
  - Pretrain epochs: 2, Finetune epochs: 2
  - Gene sampling: 500, Loyfer regions: 1125
  - Masking ratio: 0.1

**Input Data:**
- `UTERINE_TISSUE_TPM_matrix.txt` (7.9MB) - RNA expression matrix
- `UTERINE_TISSUE_meth_matrix.txt` (8.4MB) - Methylation matrix

**Output Files (MODEL_OUTPUT_V1/):**
- `best_foundation_model_*.pth` (78MB) - Best pretrained weights
- `numeric_id_prediction_model_*.pth` (78MB) - Finetuned model weights
- `metrics_*.txt` (11KB) - Training/validation metrics
- `compute_*.json` (1.2KB) - Computational statistics
- `numeric_id_predictions_*.csv` (775B) - Final predictions

**Sample Metadata (SAMPLE_IDENTIFICATION_FILES/):**
- `gdc_manifest.2025-11-11.115210.txt` (7.5KB)
- `gdc_sample_sheet.2025-11-11.tsv` (12KB)
- `metadata.repository.2025-11-11.json` (92KB)

**GDC Tools:**
- `gdc-client` (13MB binary)
- `gdc-user-token.*.txt` (547B) - Access token (⚠️ sensitive)

---

## Recommendations

### Option 1: Create Git Repository (Recommended)
```bash
cd ~/sunil_offload

# Initialize git repo
git init
git add *.py *.sh SAMPLE_IDENTIFICATION_FILES/

# Create .gitignore for large/sensitive files
cat > .gitignore <<EOF
# Large data files
*.txt
*.pth

# GDC credentials
gdc-client
gdc-user-token.*

# Output directories (regenerable)
MODEL_OUTPUT_V1/
EOF

# Commit code and configuration
git add .gitignore
git commit -m "Initial commit: RNA-methylation foundation model for TCGA deconvolution"

# Create GitHub repo and push
# (Create repo at: https://github.com/new)
git remote add origin https://github.com/MarcusAtMotley/tcga-deconvolution.git
git push -u origin main
```

### Option 2: Archive and Download
```bash
# Create archive excluding large model weights
cd ~
tar -czf sunil_offload_code.tar.gz \
  --exclude='*.pth' \
  --exclude='*.txt' \
  --exclude='gdc-client' \
  sunil_offload/

# Download to local machine
# Then optionally upload model weights separately if needed
```

### Option 3: Selective Backup
**Minimum to preserve:**
- ✅ Python scripts (code is critical)
- ✅ Shell scripts (training configuration)
- ✅ Sample metadata (GDC manifests - harder to regenerate)
- ⚠️ Model weights (156MB - regenerable but takes compute time)
- ❌ Input data matrices (16MB - available from TCGA/GDC)
- ❌ GDC client binary (13MB - downloadable)

---

## Other Directories (No Action Needed)

### ~/bin/
- **Contents:** Nextflow binary, nf-test, BaseSpace CLI (bs)
- **Action:** None needed - these are just executables, downloadable

---

## Summary

**Total Git Repositories:** 5 (all pushed ✅)
**Loose Projects Needing Backup:** 1 (sunil_offload)

**Next Steps:**
1. ✅ Decide how to backup ~/sunil_offload/ (git repo recommended)
2. ✅ Document any credentials/tokens that need to be regenerated on AWS
3. ✅ Confirm all critical code is preserved before VM deletion

**Estimated Backup Time:**
- Option 1 (Git repo creation): ~10 minutes
- Option 2 (Archive download): ~5 minutes for code, +15 min if including model weights
- Option 3 (Selective backup): ~5 minutes
