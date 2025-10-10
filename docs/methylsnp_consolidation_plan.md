# MethylSNP Pipeline Consolidation Plan

**Document Version**: 1.0
**Date**: September 24, 2025
**Author**: Development Team
**Status**: Design Phase

## Executive Summary

This document outlines a comprehensive plan to consolidate the current methylSNP analysis pipeline into a single, unified subworkflow. The primary motivation is to enable easy addition of non-hairpin processing modes while significantly simplifying the main workflow architecture.

### Key Benefits
- **Simplified main workflow**: Reduces ~60 lines of complex methylSNP orchestration to ~5 lines
- **Mode-based processing**: Easy switching between hairpin, non-hairpin, and disabled modes
- **Future extensibility**: Clean architecture for adding new processing types
- **Unified interface**: Single point of configuration and output collection
- **Improved maintainability**: Centralized methylation analysis logic

## Current Architecture Analysis

### Overview
The methylation analysis pipeline is currently distributed across multiple components:

1. **Main Workflow Orchestration** (`workflows/modulestesting.nf`)
2. **TNA_DECONVOLUTION Subworkflow** (`subworkflows/local/tna_deconvolution/main.nf`)
3. **METHYLSNP_ANALYSIS Subworkflow** (`subworkflows/local/methylsnp_analysis/main.nf`)
4. **LOFREQ_VARIANT_CALLING Subworkflow** (`subworkflows/local/lofreq_variant_calling/main.nf`)

### Current Workflow Logic (60+ lines)

```groovy
// Main workflow complexity - BEFORE
if (!params.skip_methylsnp_analysis) {
    PREPARE_REFERENCES()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    if (!params.skip_rna_deconvolution && params.rna_barcode_config) {
        // Dual processing: RNA/barcoded + DNA/unbarcoded reads

        // Process RNA/barcoded reads with transcriptome
        METHYLSNP_RNA (
            TNA_DECONVOLUTION.out.barcoded_reads,
            PREPARE_REFERENCES.out.transcriptome_index,
            PREPARE_REFERENCES.out.transcriptome_fasta
        )
        ch_versions = ch_versions.mix(METHYLSNP_RNA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_RNA.out.extraction_report.collect{it[1]})

        // Process DNA/unbarcoded reads with genome
        METHYLSNP_DNA (
            TNA_DECONVOLUTION.out.unbarcoded_reads,
            PREPARE_REFERENCES.out.genome_index,
            PREPARE_REFERENCES.out.genome_fasta
        )
        ch_versions = ch_versions.mix(METHYLSNP_DNA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_DNA.out.extraction_report.collect{it[1]})

        // Variant calling on RNA/barcoded reads
        LOFREQ_RNA (
            METHYLSNP_RNA.out.processed_bam,
            PREPARE_REFERENCES.out.transcriptome_fasta,
            PREPARE_REFERENCES.out.transcriptome_fai
        )
        ch_versions = ch_versions.mix(LOFREQ_RNA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.vcf.collect{it[1]})

        // Variant calling on DNA/unbarcoded reads
        LOFREQ_DNA (
            METHYLSNP_DNA.out.processed_bam,
            PREPARE_REFERENCES.out.genome_fasta,
            PREPARE_REFERENCES.out.genome_fai
        )
        ch_versions = ch_versions.mix(LOFREQ_DNA.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_DNA.out.vcf.collect{it[1]})

    } else {
        // Single processing: treat all reads as DNA/genome
        METHYLSNP_SINGLE (
            READ_TRIMMING.out.reads,
            PREPARE_REFERENCES.out.genome_index,
            PREPARE_REFERENCES.out.genome_fasta
        )
        ch_versions = ch_versions.mix(METHYLSNP_SINGLE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_SINGLE.out.extraction_report.collect{it[1]})

        // Variant calling on all reads (DNA/genome)
        LOFREQ_SINGLE (
            METHYLSNP_SINGLE.out.processed_bam,
            PREPARE_REFERENCES.out.genome_fasta,
            PREPARE_REFERENCES.out.genome_fai
        )
        ch_versions = ch_versions.mix(LOFREQ_SINGLE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SINGLE.out.vcf.collect{it[1]})
    }
}
```

### Current Channel Flow

```
READ_TRIMMING.out.reads
    ↓
Main Workflow Decision Point
    ↓
├─ RNA Deconvolution Mode:
│   ├─ TNA_DECONVOLUTION
│   │   ├─ barcoded_reads → METHYLSNP_RNA → LOFREQ_RNA
│   │   └─ unbarcoded_reads → METHYLSNP_DNA → LOFREQ_DNA
│   └─ MultiQC collection (scattered across 6+ mix operations)
│
└─ Single Processing Mode:
    ├─ all reads → METHYLSNP_SINGLE → LOFREQ_SINGLE
    └─ MultiQC collection (scattered across 3+ mix operations)
```

### Problems with Current Architecture

1. **Complex Main Workflow**: 60+ lines of methylSNP-specific logic clutters the main workflow
2. **Scattered MultiQC Collection**: 6+ separate `.mix()` operations spread throughout the logic
3. **Duplicated Channel Management**: Similar patterns repeated for RNA/DNA/Single paths
4. **Difficult to Extend**: Adding new processing modes requires modifying main workflow
5. **Harder to Test**: Cannot test methylation pipeline as isolated unit
6. **Version Collection Scattered**: Multiple version mix operations throughout main workflow

## Proposed Architecture: Single Consolidated Subworkflow

### New Workflow Logic (5 lines)

```groovy
// Main workflow simplicity - AFTER
if (params.methylsnp_mode != 'disabled') {
    COMPREHENSIVE_METHYLSNP_ANALYSIS(
        READ_TRIMMING.out.reads,
        PREPARE_REFERENCES.out.genome_fasta,
        PREPARE_REFERENCES.out.transcriptome_fasta,
        PREPARE_REFERENCES.out.genome_index,
        PREPARE_REFERENCES.out.transcriptome_index,
        PREPARE_REFERENCES.out.genome_fai,
        PREPARE_REFERENCES.out.transcriptome_fai,
        params.methylsnp_mode,  // 'hairpin', 'non_hairpin', 'disabled'
        params.rna_barcode_config
    )
    ch_versions = ch_versions.mix(COMPREHENSIVE_METHYLSNP_ANALYSIS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COMPREHENSIVE_METHYLSNP_ANALYSIS.out.multiqc_files)
}
```

### Internal Subworkflow Structure

```groovy
workflow COMPREHENSIVE_METHYLSNP_ANALYSIS {

    take:
    ch_reads              // channel: [meta, [fastq1, fastq2]]
    genome_fasta          // channel: [meta, genome.fa]
    transcriptome_fasta   // channel: [meta, transcriptome.fa]
    genome_index          // channel: [meta, [*.bt2]]
    transcriptome_index   // channel: [meta, [*.bt2]]
    genome_fai            // channel: [meta, genome.fa.fai]
    transcriptome_fai     // channel: [meta, transcriptome.fa.fai]
    mode                  // val: 'hairpin', 'non_hairpin'
    rna_config            // path: RNA barcode config file (optional)

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (mode == 'hairpin') {
        //
        // HAIRPIN MODE: Current pipeline logic
        //
        if (rna_config && !params.skip_rna_deconvolution) {
            // RNA/DNA dual processing
            ch_reads_with_cfg = ch_reads.map { meta, reads ->
                tuple(meta, reads, file(rna_config))
            }

            TNA_DECONVOLUTION(ch_reads_with_cfg)
            ch_versions = ch_versions.mix(TNA_DECONVOLUTION.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.txt_report.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.json_report.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.rna_stats_mqc)

            // RNA/barcoded processing
            METHYLSNP_ANALYSIS as METHYLSNP_RNA (
                TNA_DECONVOLUTION.out.barcoded_reads,
                transcriptome_index,
                transcriptome_fasta
            )
            ch_versions = ch_versions.mix(METHYLSNP_RNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_RNA.out.extraction_report.collect{it[1]})

            LOFREQ_VARIANT_CALLING as LOFREQ_RNA (
                METHYLSNP_RNA.out.processed_bam,
                transcriptome_fasta,
                transcriptome_fai
            )
            ch_versions = ch_versions.mix(LOFREQ_RNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.vcf.collect{it[1]})

            // DNA/unbarcoded processing
            METHYLSNP_ANALYSIS as METHYLSNP_DNA (
                TNA_DECONVOLUTION.out.unbarcoded_reads,
                genome_index,
                genome_fasta
            )
            ch_versions = ch_versions.mix(METHYLSNP_DNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_DNA.out.extraction_report.collect{it[1]})

            LOFREQ_VARIANT_CALLING as LOFREQ_DNA (
                METHYLSNP_DNA.out.processed_bam,
                genome_fasta,
                genome_fai
            )
            ch_versions = ch_versions.mix(LOFREQ_DNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_DNA.out.vcf.collect{it[1]})

            // Combine outputs for unified interface
            ch_processed_bam = METHYLSNP_RNA.out.processed_bam.mix(METHYLSNP_DNA.out.processed_bam)
            ch_vcf_files = LOFREQ_RNA.out.vcf.mix(LOFREQ_DNA.out.vcf)

        } else {
            // Single processing mode
            METHYLSNP_ANALYSIS as METHYLSNP_SINGLE (
                ch_reads,
                genome_index,
                genome_fasta
            )
            ch_versions = ch_versions.mix(METHYLSNP_SINGLE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_SINGLE.out.extraction_report.collect{it[1]})

            LOFREQ_VARIANT_CALLING as LOFREQ_SINGLE (
                METHYLSNP_SINGLE.out.processed_bam,
                genome_fasta,
                genome_fai
            )
            ch_versions = ch_versions.mix(LOFREQ_SINGLE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SINGLE.out.vcf.collect{it[1]})

            // Single outputs
            ch_processed_bam = METHYLSNP_SINGLE.out.processed_bam
            ch_vcf_files = LOFREQ_SINGLE.out.vcf
        }

    } else if (mode == 'non_hairpin') {
        //
        // NON-HAIRPIN MODE: Future implementation
        //
        // TODO: Implement alternative processing pipeline
        // - Skip TNA deconvolution entirely
        // - Use different trimming/processing parameters
        // - Alternative alignment strategy
        // - Same variant calling endpoint

        // Placeholder for future implementation:
        // NON_HAIRPIN_METHYLSNP_ANALYSIS(ch_reads, genome_index, genome_fasta)
        // LOFREQ_VARIANT_CALLING(...)

        error "Non-hairpin mode not yet implemented"
    }

    emit:
    processed_bam  = ch_processed_bam   // channel: [meta, *.bam]
    vcf_files      = ch_vcf_files       // channel: [meta, *.vcf.gz]
    multiqc_files  = ch_multiqc_files   // channel: [*]
    versions       = ch_versions        // channel: [versions.yml]
}
```

### New Channel Flow

```
READ_TRIMMING.out.reads
    ↓
COMPREHENSIVE_METHYLSNP_ANALYSIS
    ├─ Mode Selection (hairpin/non_hairpin)
    ├─ Internal processing (all current logic)
    ├─ Unified MultiQC collection
    └─ Unified version collection
    ↓
Clean unified outputs
```

## Parameter Design Enhancement

### Current Parameters
```groovy
--skip_methylsnp_analysis     // Binary: true/false
--skip_rna_deconvolution      // Binary: true/false
--rna_barcode_config          // Path to config file
```

### Proposed Parameters
```groovy
--methylsnp_mode              // String: 'hairpin', 'non_hairpin', 'disabled'
--rna_barcode_config          // Path to config file (for hairpin mode)
--non_hairpin_config          // Path to config file (for non-hairpin mode)
```

### Parameter Logic
- `methylsnp_mode = 'disabled'` → Skip all methylation analysis
- `methylsnp_mode = 'hairpin'` → Current pipeline behavior
- `methylsnp_mode = 'non_hairpin'` → Future alternative pipeline

## Implementation Roadmap

### Phase 1: Consolidation (Current Functionality)
**Estimated Time**: 2-3 days

1. **Create new subworkflow structure**
   - `subworkflows/local/comprehensive_methylsnp_analysis/main.nf`
   - Move all current hairpin logic from main workflow

2. **Update main workflow**
   - Replace 60+ lines with single subworkflow call
   - Preserve all current functionality

3. **Add new parameter**
   - Add `methylsnp_mode` to `nextflow_schema.json`
   - Default to `'hairpin'` for backward compatibility
   - Map old `skip_methylsnp_analysis` to `'disabled'` mode

4. **Testing and validation**
   - Verify identical behavior to current pipeline
   - Test all three modes: RNA deconvolution, single processing, disabled

### Phase 2: Non-Hairpin Mode Framework (2-4 weeks)
**Estimated Time**: 1-2 weeks

1. **Design non-hairpin processing requirements**
   - Define what differs from hairpin processing
   - Identify which modules need alternatives
   - Plan parameter requirements

2. **Implement non-hairpin branch**
   - Add `mode == 'non_hairpin'` branch to subworkflow
   - Create placeholder implementations
   - Add proper error handling for unimplemented features

3. **Testing framework**
   - Add mode-specific test cases
   - Validation for each mode independently

### Phase 3: Non-Hairpin Implementation (Timeline TBD)
**Estimated Time**: Variable based on requirements

1. **Implement non-hairpin specific modules**
2. **Add non-hairpin specific parameters**
3. **Full testing and validation**
4. **Documentation updates**

## Benefits Analysis

### Immediate Benefits (Phase 1)
- **Reduced Main Workflow Complexity**: 60+ lines → 5 lines
- **Unified MultiQC Collection**: Single point instead of 6+ scattered operations
- **Unified Version Collection**: Single point instead of multiple scattered operations
- **Improved Testability**: Can test methylation pipeline in isolation
- **Better Error Handling**: Failures contained within subworkflow
- **Cleaner Abstractions**: Main workflow focuses on high-level orchestration

### Future Benefits (Phase 2-3)
- **Easy Mode Addition**: Simple to add new processing types
- **Mode-Specific Parameters**: Clean parameter validation per mode
- **Parallel Development**: Different modes can be developed independently
- **Mixed Processing**: Future possibility of per-sample mode selection

## Risk Assessment and Mitigation

### Implementation Risks

**Risk**: Breaking existing functionality during consolidation
- **Mitigation**: Implement alongside current structure for parallel testing
- **Rollback**: Keep original structure until new version fully validated

**Risk**: Complex channel management within subworkflow
- **Mitigation**: Extensive testing of channel flows
- **Validation**: Use same test patterns as current pipeline

**Risk**: Performance impact from additional subworkflow layer
- **Mitigation**: Profile before/after performance
- **Expectation**: Minimal impact, possibly improved due to better channel management

### Compatibility Risks

**Risk**: Parameter changes breaking existing configurations
- **Mitigation**: Maintain backward compatibility with parameter mapping
- **Transition**: Gradual deprecation of old parameters with warnings

**Risk**: Changes to output structure
- **Mitigation**: Maintain identical output channels and file structures
- **Validation**: Comprehensive output testing

## Testing Strategy

### Unit Testing
- **Subworkflow Testing**: Test COMPREHENSIVE_METHYLSNP_ANALYSIS independently
- **Mode Testing**: Test each mode (hairpin/non_hairpin/disabled) separately
- **Channel Testing**: Validate all input/output channel structures

### Integration Testing
- **Main Workflow**: Test full pipeline with new subworkflow
- **Comparison Testing**: Compare outputs between old and new architecture
- **Edge Case Testing**: Test with various parameter combinations

### Regression Testing
- **Existing Datasets**: Run on known good datasets
- **Output Validation**: Ensure identical results to current pipeline
- **Performance Testing**: Validate no significant performance degradation

## Technical Considerations

### Channel Management
- **Input Channels**: Clean parameter passing to subworkflow
- **Internal Channels**: Proper channel scoping within subworkflow
- **Output Channels**: Unified emit structure regardless of internal branching

### Error Handling
- **Mode Validation**: Proper validation of mode parameter
- **Graceful Failures**: Clean error messages for configuration issues
- **Partial Failures**: Handle failures in one branch without affecting others

### Resource Management
- **Memory Usage**: Monitor memory usage patterns in consolidated subworkflow
- **CPU Usage**: Ensure efficient resource utilization
- **Disk Usage**: Consider disk space for intermediate files

### Extensibility
- **New Modes**: Framework for adding additional processing modes
- **Parameter Expansion**: Structure for mode-specific parameters
- **Module Integration**: Easy integration of new processing modules

## Migration Strategy

### Backward Compatibility
1. **Parameter Mapping**: Map old parameters to new mode system
2. **Deprecation Warnings**: Warn users about parameter changes
3. **Documentation**: Update all documentation with new parameter usage

### Rollback Plan
1. **Parallel Implementation**: Keep both architectures during transition
2. **Feature Flag**: Allow switching between old/new architecture
3. **Validation**: Extensive testing before removing old architecture

## Future Considerations

### Non-Hairpin Processing Requirements
- **Alternative Trimming**: Different adapter trimming strategies
- **Alternative Alignment**: Different alignment parameters or tools
- **Alternative Processing**: Skip deconvolution-specific steps
- **Alternative QC**: Different quality control metrics

### Advanced Features
- **Per-Sample Modes**: Allow different samples to use different modes
- **Mixed Processing**: Process some samples as hairpin, others as non-hairpin
- **Dynamic Mode Selection**: Auto-detect appropriate mode based on sample characteristics

### Monitoring and Observability
- **Mode Tracking**: Track which mode is being used in logs
- **Performance Metrics**: Monitor performance by mode
- **Usage Analytics**: Track adoption of different modes

## Conclusion

The consolidation of methylSNP processing into a single subworkflow represents a significant architectural improvement that:

1. **Simplifies the main workflow** from 60+ lines to 5 lines
2. **Enables easy addition of non-hairpin processing** through clean mode switching
3. **Improves maintainability** by centralizing all methylation logic
4. **Enhances testability** by creating isolated, testable components
5. **Provides a foundation for future extensibility** with additional processing modes

This approach directly addresses the goal of making methylSNP steps "easily skippable" while creating a clean path for adding non-hairpin processing alternatives. The phased implementation approach ensures minimal risk while maximizing the architectural benefits.

---

**Next Steps**: Review this plan with the development team and proceed with Phase 1 implementation upon approval.