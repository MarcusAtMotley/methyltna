# methyltna Pipeline - Future Enhancements & Module Roadmap

**Document Version**: 1.0
**Last Updated**: October 10, 2025
**Current Pipeline Version**: 1.0.0dev

---

## 🎯 Executive Summary

This document outlines potential enhancements to the methyltna pipeline. The recommendations are organized by scientific impact and implementation priority, with special emphasis on leveraging the unique advantages of TNA-EM-seq technology.

### Top 3 Strategic Priorities

1. **Multi-Omics Integration** - Correlate RNA expression with DNA methylation (unique to TNA-EM-seq!)
2. **Differential Analysis Suite** - Statistical analysis for both RNA and methylation data
3. **Variant Annotation** - Functional interpretation of detected variants

---

## 🔬 High-Priority Scientific Modules

### 1. Differential Expression Analysis (RNA-seq)

**Priority**: ⭐⭐⭐⭐⭐
**Implementation Effort**: Medium
**Scientific Impact**: High

#### Description
Add statistical analysis of gene expression data from FeatureCounts output.

#### Recommended Tools
- **DESeq2** (R-based, most popular)
- **edgeR** (R-based, flexible)
- **limma-voom** (R-based, handles weights well)

#### Implementation Plan
- Create R-based nf-core module or custom module
- Input: FeatureCounts TSV files
- Requires: Experimental design/comparison groups
- Output formats:
  - DE gene lists (CSV/TSV)
  - MA plots
  - Volcano plots
  - PCA plots
  - Heatmaps

#### MultiQC Integration
- DE summary statistics table
- Top DE genes visualization
- QC metrics (dispersion, normalization factors)

#### Benefits
- Complete RNA-seq analysis workflow
- Publication-ready figures and statistics
- Biological interpretation of expression changes

---

### 2. Differential Methylation Analysis

**Priority**: ⭐⭐⭐⭐⭐
**Implementation Effort**: Medium-High
**Scientific Impact**: High

#### Description
Identify differentially methylated regions (DMRs) between experimental conditions.

#### Recommended Tools
- **DSS** (R-based, handles biological replicates well)
- **methylKit** (R-based, comprehensive)
- **dmrseq** (R-based, region-based approach)

#### Implementation Plan
- R-based module processing Biscuit VCF outputs
- Input: Methylation calls from BISCUIT_PILEUP
- Requires: Sample grouping/comparisons
- Output formats:
  - DMR lists (BED files)
  - Methylation heatmaps
  - Circos plots (genome-wide view)
  - Methylation profiles around genes/TSSs

#### MultiQC Integration
- DMR summary statistics
- Methylation distribution plots
- Context-specific methylation (CpG, CHG, CHH)

#### Benefits
- Complete methylation analysis workflow
- Identify epigenetically regulated regions
- Link methylation changes to biological phenotypes

---

### 3. Multi-Omics Integration ⭐ UNIQUE VALUE!

**Priority**: ⭐⭐⭐⭐⭐
**Implementation Effort**: High
**Scientific Impact**: VERY HIGH

#### Description
**This is the key advantage of TNA-EM-seq** - integrate RNA expression with DNA methylation from the same sample.

#### Analysis Capabilities
1. **Promoter Methylation vs Expression**
   - Correlate gene expression with upstream methylation
   - Identify methylation-silenced genes
   - Detect inverse correlation patterns

2. **Allele-Specific Analysis**
   - Link methylation patterns to allelic expression
   - Identify imprinted genes
   - Detect allele-specific methylation

3. **Regulatory Element Analysis**
   - Correlate enhancer methylation with target gene expression
   - Identify tissue-specific regulatory mechanisms
   - Map methylation-expression networks

#### Recommended Tools
- **Custom R/Python pipeline** integrating:
  - FeatureCounts output (RNA)
  - Biscuit methylation calls (DNA)
  - LoFreq variants (for phasing)
- **MOFA2** (Multi-Omics Factor Analysis)
- **mixOmics** (R package for multi-omics integration)

#### Output Formats
- Correlation matrices (gene-wise, region-wise)
- Integrated heatmaps (expression + methylation)
- Network graphs (methylation-expression relationships)
- Pathway enrichment (genes with coordinated changes)
- Interactive HTML reports

#### Benefits
- **Unlock unique TNA-EM-seq advantages**
- Identify epigenetic regulation mechanisms
- Discover novel biomarkers (expression+methylation signatures)
- Understand cause-effect relationships (methylation → expression)

---

### 4. Variant Annotation

**Priority**: ⭐⭐⭐⭐
**Implementation Effort**: Low-Medium
**Scientific Impact**: Medium-High

#### Description
Annotate functional impact of variants detected by LoFreq.

#### Recommended Tools
- **VEP (Variant Effect Predictor)** - Ensembl
- **SnpEff** - Faster, self-contained
- **ANNOVAR** - Comprehensive annotations

#### Implementation Plan
- nf-core module: `ensemblvep` or custom SnpEff module
- Input: LoFreq VCF files
- Output formats:
  - Annotated VCFs
  - Impact summaries (HIGH/MODERATE/LOW)
  - Gene-level variant reports
  - Pathogenicity predictions (SIFT, PolyPhen)

#### MultiQC Integration
- Variant impact distribution plots
- Affected gene lists
- Variant type summaries

#### Benefits
- Understand biological significance of variants
- Prioritize variants for validation
- Identify disease-relevant mutations
- Link variants to expression/methylation changes

---

## 🔧 Technical Improvements

### 5. Replace MarkDup.py with Picard MarkDuplicates

**Priority**: ⭐⭐⭐
**Implementation Effort**: Low
**Scientific Impact**: Medium (maintainability)

#### Current Implementation
- Custom Python 2 script in METHYLSNP_PROCESSING
- Groups by: chromosome + position + sequence identity
- Keeps highest MAPQ read

#### Proposed Implementation
- Use nf-core `picard/markduplicates` module
- Industry standard, well-tested
- Much faster (Java vs Python 2)
- Already generates MultiQC-compatible metrics

#### Key Differences
- **Picard grouping**: Library + Reference + Strand + 5' position (NO sequence identity check)
- **Picard selection**: Sum of base qualities (not MAPQ)
- **Potential impact**: May be more aggressive in marking duplicates

#### Testing Required
1. Install module: `nf-core modules install picard/markduplicates`
2. Run comparison: MarkDup.py vs Picard on test dataset
3. Compare read counts and methylation calling results
4. **Validate**: Is sequence-identity check critical for methylSNP?

#### Benefits
- Python 2 deprecation resolved
- Faster processing
- Better QC metrics
- Industry-standard approach

---

### 6. UMI (Unique Molecular Identifier) Support

**Priority**: ⭐⭐⭐
**Implementation Effort**: Medium
**Scientific Impact**: High (if UMIs used in protocol)

#### Description
If your sequencing protocol includes UMIs, use them for accurate deduplication.

#### Recommended Tools
- **UMI-tools** (Python-based, widely used)
- **fgbio** (Scala-based, Broad Institute)

#### Implementation Plan
- Add UMI extraction step (before/after trimming)
- UMI-aware deduplication (replaces MarkDuplicates)
- UMI consensus calling for error correction

#### Benefits
- More accurate quantification (true biological vs PCR duplicates)
- Better variant calling (consensus reduces errors)
- Essential for low-input samples

---

### 7. Allele-Specific Expression & Methylation Analysis

**Priority**: ⭐⭐⭐⭐
**Implementation Effort**: High
**Scientific Impact**: Very High

#### Description
Detect allele-specific patterns in both RNA and DNA using phased variants.

#### Recommended Tools
- **WASP** (allele-specific alignment filtering)
- **phASER** (phasing and ASE analysis)
- **Custom pipeline** leveraging LoFreq variants

#### Implementation Plan
1. Phase variants from LoFreq output
2. Separate reads by allele (WASP filtering)
3. Quantify allele-specific expression (RNA)
4. Quantify allele-specific methylation (DNA)
5. Correlate ASE with allelic methylation

#### Outputs
- ASE gene lists
- Allelic imbalance statistics
- Methylation differences by allele
- Imprinted gene identification

#### Benefits
- **Perfect fit for TNA-EM-seq!** (variants + RNA + methylation)
- Identify cis-regulatory variants
- Detect imprinting and X-inactivation
- Understand allelic regulation

---

## 📊 Quality & Reporting Enhancements

### 8. Interactive Visualization

**Priority**: ⭐⭐⭐
**Implementation Effort**: Medium
**Scientific Impact**: Medium (usability)

#### Description
Browser-based visualization of alignments and methylation patterns.

#### Recommended Tools
- **IGV.js** (JavaScript IGV embedded in HTML)
- **JBrowse 2** (next-gen genome browser)
- **R Shiny** (interactive R-based apps)

#### Implementation Plan
- Generate IGV.js HTML reports
- Create JBrowse track hub
- Embed in MultiQC or separate reports

#### Outputs
- Standalone HTML genome browser
- Pre-loaded tracks (alignment, methylation, variants)
- Interactive exploration (zoom, pan, search)

#### Benefits
- No IGV installation required
- Easy sharing with collaborators
- Publication-quality screenshots

---

### 9. Contamination Screening

**Priority**: ⭐⭐⭐
**Implementation Effort**: Low
**Scientific Impact**: Medium (QC)

#### Description
Detect contaminating sequences from other organisms.

#### Recommended Tools
- **FastQ Screen** (multi-genome screening)
- **Kraken2** (taxonomic classification)
- **Centrifuge** (metagenomic classification)

#### Implementation Plan
- Add screening step after read trimming
- Screen against common contaminants (bacteria, fungi, viral, etc.)
- Report contamination percentages

#### MultiQC Integration
- Contamination summary plots
- Species composition charts
- Flagging high-contamination samples

#### Benefits
- Quality assurance
- Identify experimental contamination
- Validate sample purity

---

### 10. Strand-Specific Metrics

**Priority**: ⭐⭐
**Implementation Effort**: Low
**Scientific Impact**: Medium (QC)

#### Description
Validate strand-specificity of RNA-seq library preparation.

#### Recommended Tools
- **RSeQC** (comprehensive RNA-seq QC)
- **Qualimap** (alignment QC with strand metrics)

#### Implementation Plan
- Add RSeQC `infer_experiment.py`
- Check strand balance in STAR alignments
- Generate junction saturation curves

#### MultiQC Integration
- Strand balance plots
- Junction saturation curves
- Read distribution (exon/intron/intergenic)

#### Benefits
- Validate library prep protocol
- Detect strand-swapping issues
- Ensure proper directionality

---

## 🧬 Advanced Analysis Modules

### 11. Structural Variant Detection

**Priority**: ⭐⭐
**Implementation Effort**: Medium-High
**Scientific Impact**: High (for cancer/disease studies)

#### Description
Detect large-scale genomic rearrangements (deletions, duplications, inversions, translocations).

#### Recommended Tools
- **Manta** (Illumina, fast)
- **GRIDSS** (comprehensive)
- **Delly** (balanced SV detection)

#### Implementation Plan
- Use Biscuit DNA alignments as input
- Call SVs with paired-end/split-read evidence
- Filter and annotate SV breakpoints

#### Outputs
- SV VCF files
- Circos plots (genome-wide view)
- Gene disruption analysis

#### Benefits
- Comprehensive variant landscape
- Identify cancer-relevant rearrangements
- Detect gene fusions

---

### 12. Copy Number Variation (CNV) Analysis

**Priority**: ⭐⭐
**Implementation Effort**: Medium
**Scientific Impact**: High (for cancer studies)

#### Description
Detect chromosomal gains and losses from read depth.

#### Recommended Tools
- **Control-FREEC** (no control samples needed)
- **CNVkit** (targeted/WGS)
- **ichorCNA** (cell-free DNA)

#### Implementation Plan
- Use Biscuit DNA alignment read depth
- Normalize for GC bias and mappability
- Segment genome and call CNVs

#### Outputs
- CNV segments (gains/losses)
- Genome-wide CNV plots
- Gene-level copy number

#### Benefits
- Identify genomic instability
- Detect oncogene amplifications
- Tumor suppressor deletions

---

### 13. Splicing Analysis

**Priority**: ⭐⭐⭐
**Implementation Effort**: Medium
**Scientific Impact**: High

#### Description
Detect differential splicing events and alternative isoforms.

#### Recommended Tools
- **rMATS** (replicate MATS, widely used)
- **LeafCutter** (intron-centric approach)
- **SUPPA2** (fast, flexible)

#### Implementation Plan
- Use STAR alignments (already generated!)
- Compare splicing between conditions
- Classify event types (exon skip, alt 5'/3'SS, etc.)

#### Outputs
- Differential splicing events
- Sashimi plots (junction visualization)
- Isoform switching analysis

#### Benefits
- **Leverage existing STAR alignments**
- Understand post-transcriptional regulation
- Identify disease-relevant splicing changes

---

## 🚀 Workflow Enhancements

### 14. Batch Effect Correction

**Priority**: ⭐⭐⭐
**Implementation Effort**: Low-Medium
**Scientific Impact**: High (multi-batch studies)

#### Description
Correct technical variation across sequencing runs or batches.

#### Recommended Tools
- **ComBat-seq** (count-based correction)
- **RUVseq** (Remove Unwanted Variation)
- **sva** (Surrogate Variable Analysis)

#### Implementation Plan
- Apply before differential analysis
- Requires batch information in samplesheet
- Preserve biological variation

#### Benefits
- More robust statistical results
- Essential for multi-center studies
- Reduce false positives

---

### 15. Pathway Enrichment Analysis

**Priority**: ⭐⭐⭐⭐
**Implementation Effort**: Low-Medium
**Scientific Impact**: High (interpretation)

#### Description
Biological interpretation of DE genes and DMRs.

#### Recommended Tools
- **clusterProfiler** (R, comprehensive)
- **GSEA** (Gene Set Enrichment Analysis)
- **Enrichr** (web-based, many databases)
- **g:Profiler** (multi-species support)

#### Implementation Plan
- Input: DE gene lists, DMR-associated genes
- Query: GO, KEGG, Reactome, MSigDB
- Generate: Enrichment plots, networks

#### Outputs
- Enriched pathway tables
- Dotplots, barplots
- Enrichment maps
- Gene-concept networks

#### Benefits
- Understand biological mechanisms
- Identify dysregulated pathways
- Generate hypotheses for follow-up

---

## 📋 Implementation Priority Matrix

### Immediate Priorities (Next 3-6 months)
1. ✅ **Multi-Omics Integration** - Core value of TNA-EM-seq
2. ✅ **Differential Analysis** (DESeq2 + DSS/methylKit)
3. ✅ **Variant Annotation** (VEP/SnpEff)
4. ⚙️ **Pathway Enrichment** (clusterProfiler)

### Medium-Term (6-12 months)
5. ⚙️ **Allele-Specific Analysis** (WASP/phASER)
6. ⚙️ **Replace MarkDup.py** with Picard
7. ⚙️ **Splicing Analysis** (rMATS/LeafCutter)
8. ⚙️ **Interactive Visualization** (IGV.js)

### Long-Term (12+ months)
9. 📅 **Structural Variants** (Manta/GRIDSS)
10. 📅 **CNV Analysis** (Control-FREEC)
11. 📅 **UMI Support** (if protocol updated)
12. 📅 **Contamination Screening** (Kraken2)

---

## 🛠️ Technical Considerations

### Container Support
- Ensure all new modules have:
  - Conda environments
  - Docker containers
  - Singularity/Apptainer support

### Resource Requirements
- **R-based modules**: 16-32GB RAM for DE analysis
- **Annotation modules**: 8-16GB RAM + reference databases
- **SV/CNV callers**: 32-64GB RAM for human genome

### Testing Strategy
- Create nf-test for each new module
- Add to CI/CD pipeline
- Validate with test datasets
- Compare results to standalone tool runs

### Documentation Needs
- Update README with new features
- Add usage examples
- Create tutorial notebooks
- Update pipeline flowchart

---

## 💡 Research-Specific Recommendations

### For Cancer Studies
**Priority modules:**
1. CNV Analysis (genomic instability)
2. Structural Variants (gene fusions)
3. Multi-Omics Integration (methylation-silencing)
4. Variant Annotation (driver mutations)

### For Developmental Biology
**Priority modules:**
1. Differential Methylation (cell type differences)
2. Allele-Specific Analysis (imprinting)
3. Multi-Omics Integration (epigenetic regulation)
4. Splicing Analysis (isoform switching)

### For Disease Biomarker Discovery
**Priority modules:**
1. Multi-Omics Integration (combined signatures)
2. Differential Analysis (both RNA + methylation)
3. Pathway Enrichment (mechanism understanding)
4. Interactive Visualization (presentation)

---

## 📚 Resources & References

### Useful Links
- [nf-core modules](https://nf-co.re/modules) - Pre-built Nextflow modules
- [Bioconductor](https://bioconductor.org/) - R packages for analysis
- [MultiQC plugin docs](https://multiqc.info/docs/) - Custom reporting

### Key Papers
- TNA-EM-seq methodology (if published)
- Multi-omics integration approaches
- Methylation analysis best practices
- RNA-seq analysis guidelines

### Training Resources
- nf-core/bytesize talks
- Bioconductor workflows
- Galaxy training materials

---

## 🤝 Community Contributions

### How to Contribute
1. Pick a module from this roadmap
2. Check nf-core modules repository for existing solutions
3. Follow nf-core module development guidelines
4. Submit PR with tests and documentation

### Contact & Discussion
- GitHub Issues: Feature requests and bugs
- Discussions: Design decisions and ideas
- Slack/Discord: Real-time collaboration

---

*This is a living document. Please update as new priorities emerge and modules are implemented.*

**Next Review Date**: January 2026
