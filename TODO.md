# TODO - RNA Barcode Extraction Integration

## ✅ Completed Tasks

### RNA Barcode Extraction Module (September 15, 2025)
- [x] **Created RNABARCODEEXTRACTION module** in `modules/local/rnabarcodeextraction/`
- [x] **Integrated portable cutadapt wrapper** (`portable_cutadapt_wrapper.py`) in templates directory
- [x] **Fixed version checking** - handles `--version` flag properly
- [x] **Fixed output file patterns** - matches actual cutadapt output structure
- [x] **Tested container compatibility**:
  - [x] Conda environment with `uv>=0.4.0`
  - [x] Apptainer/Singularity with `ghcr.io/astral-sh/uv:python3.12-bookworm`
- [x] **Verified functionality** with test data:
  - Sample: Ca3-RNA-NNSR_09Z_S33_L001 (501,083 read pairs)
  - Results: 55.1% barcoded, 19.4% double-barcoded, 25.5% unbarcoded
  - Outputs: barcoded/unbarcoded FASTQ files, JSON/text reports

## 🎯 Pending Tasks

### High Priority - Subworkflow Integration
- [ ] **Create barcode extraction subworkflow** using `nf-core subworkflows create`
- [ ] **Integrate barcode extraction** into main pipeline workflow
- [ ] **Add barcode extraction parameters** to `nextflow.config`
  - [ ] `skip_barcode_extraction` (boolean, default: false)
  - [ ] `barcode_config_file` (path to YAML config)
  - [ ] `barcode_adapter_name` (string, default: 'NNSR')
  - [ ] `barcode_adapter_sequence` (string)
  - [ ] `barcode_error_rate` (float, default: 0.0)
  - [ ] `barcode_min_overlap` (int, default: 10)
  - [ ] `barcode_action` (string, default: 'lowercase')
- [ ] **Update `nextflow_schema.json`** with parameter definitions and validation
- [ ] **Update modules.config** with barcode extraction process settings

### Medium Priority - Workflow Enhancement
- [ ] **See if MultiQC has a setting** to keep the different FastQC runs in separate tabs in the report
- [ ] **Add barcode extraction** to MultiQC reporting
- [ ] **Create test profile** for barcode extraction
- [ ] **Update samplesheet schema** if needed for barcode-specific metadata
- [ ] **Add barcode extraction** to pipeline documentation

### Low Priority - Advanced Features
- [ ] **Create barcode extraction tests** using nf-test
- [ ] **Add barcode QC plots** to MultiQC
- [ ] **Support multiple barcode configs** for different sample types
- [ ] **Add barcode statistics** to pipeline summary

## 📁 File Structure Created

```
modules/local/rnabarcodeextraction/
├── main.nf                           # ✅ Module definition
├── meta.yml                          # ✅ Module metadata
├── environment.yml                   # ✅ Conda dependencies (uv only)
├── templates/
│   └── portable_cutadapt_wrapper.py  # ✅ Cutadapt wrapper script
└── tests/
    └── main.nf.test                  # ✅ Basic test structure
```

## 🧪 Test Files Available

- **Config**: `runner_settings_readStart_EM_e0.0.yaml` (NNSR barcode extraction settings)
- **Test script**: `test_barcode_module.nf` (standalone module test)
- **Test data**: RNA FASTQ files in `/home/marcus/projects/mot23/nf_output/trimgalore/`

## 📝 Next Session Goals

1. **Create subworkflow** using nf-core tools
2. **Add parameters** to pipeline configuration
3. **Test integration** with existing READ_TRIMMING workflow
4. **Verify MultiQC** includes barcode extraction reports

## 🔧 Commands for Next Session

```bash
# Create subworkflow
nf-core subworkflows create barcode_extraction

# Test integrated pipeline
nextflow run . -profile test,conda --skip_barcode_extraction false --outdir results

# Validate configuration
nextflow config -profile test
```

---
*Last updated: September 15, 2025*
*Module ready for integration - all core functionality tested and working*