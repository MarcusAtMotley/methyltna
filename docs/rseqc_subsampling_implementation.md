# RSeQC Subsampling Implementation Plan

## Background

RSeQC quality control analyses, particularly `geneBody_coverage.py`, are very time-consuming when run on full-depth BAM files. For a 20-sample mixed run with ~254K transcripts:
- Current runtime: **~40 hours total** with extreme variability by sample type:
  - DNA samples (4M reads): ~24 min/sample (~5.7 min per million reads)
  - TNA samples (35M reads): ~183 min/sample (~5.4 min per million reads)
  - RNA samples (10M reads): ~86 min/sample (~11.1 min per million reads)
- **RNA samples take 2x longer per million reads** than DNA/TNA due to better transcript coverage
- Most time spent on geneBody_coverage step
- Other RSeQC tools (bamstat, read_distribution, inner_distance) are relatively fast

## Problem Statement

The pipeline currently processes heterogeneous samples with varying sequencing depths:
- High coverage samples: 30M+ reads → very slow QC
- Low coverage samples: <1M reads → acceptable speed but potentially noisy metrics
- **Challenge**: Need consistent, comparable QC metrics across all samples regardless of depth

## Solution: Fixed Read Count Subsampling

Subsample BAM files to a **fixed maximum read count** before running RSeQC analyses.

### Why Fixed Count vs. Percentage?

**Percentage-based (NOT recommended):**
- 30M reads × 10% = 3M reads ✓
- 1M reads × 10% = 100K reads ✗ (too few, noisy)
- Inconsistent QC depth across samples

**Fixed count (RECOMMENDED):**
- 30M reads → subsample to 5M reads (6x speedup)
- 1M reads → use all 1M reads (no loss)
- Consistent, comparable QC metrics across all samples

## Technical Specifications

### Subsampling Parameters

**Recommended threshold: 5 million mapped reads**
- Provides excellent statistical power for all RSeQC metrics
- Caps runtime while maintaining quality
- Still allows fair comparison of lower-depth samples

Alternative thresholds:
- **2M reads**: Faster, still reliable for most QC metrics
- **10M reads**: More conservative, but slower

### Implementation Details

1. **Count mapped reads** (exclude unmapped and secondary alignments):
   ```bash
   samtools view -c -F 260 input.bam
   ```
   - Flag 260 = unmapped (4) + not primary alignment (256)

2. **Calculate subsampling fraction**:
   ```bash
   fraction = max_reads / total_mapped_reads
   ```

3. **Subsample using samtools**:
   ```bash
   samtools view -s $fraction -b input.bam > subsampled.bam
   ```
   - The `-s` flag uses seed + fraction (e.g., `-s 0.42` = 42% random sample)
   - Random sampling ensures representative coverage

4. **Run RSeQC on subsampled BAMs**

5. **Clean up temporary subsampled BAMs** (optional)

## Code Implementation

### Bash Function

```bash
#!/bin/bash
# Subsample BAM file to maximum read count for QC
# Usage: subsample_bam_for_qc <input.bam> <output.bam> [max_reads]

subsample_bam_for_qc() {
    local input_bam=$1
    local output_bam=$2
    local max_reads=${3:-5000000}  # Default: 5M reads

    # Count mapped reads (exclude unmapped and secondary)
    local total_reads=$(samtools view -c -F 260 "${input_bam}")

    echo "Sample: $(basename ${input_bam})"
    echo "  Total mapped reads: ${total_reads}"

    if [ $total_reads -gt $max_reads ]; then
        # Calculate subsampling fraction
        local fraction=$(awk "BEGIN {printf \"%.6f\", ${max_reads}/${total_reads}}")
        echo "  Subsampling to ${max_reads} reads (${fraction}x)"

        # Subsample
        samtools view -s ${fraction} -b "${input_bam}" > "${output_bam}"

        # Index subsampled BAM
        samtools index "${output_bam}"
    else
        # Use full BAM (below threshold)
        echo "  Using all reads (below ${max_reads} threshold)"
        cp "${input_bam}" "${output_bam}"
        cp "${input_bam}.bai" "${output_bam}.bai" 2>/dev/null || samtools index "${output_bam}"
    fi
}

# Example usage:
# subsample_bam_for_qc input.sorted.bam temp.subsampled.bam 5000000
```

### Integration with RSeQC Script

Modified workflow:
1. For each BAM file
2. Create subsampled version (if needed)
3. Run RSeQC on subsampled BAM
4. Clean up subsampled BAM
5. Move to next sample

```bash
for bam in results/alignment/sorted_bam/*.bam; do
    sample=$(basename ${bam} .sorted.bam)

    # Create temporary subsampled BAM
    temp_bam="temp_${sample}_subsampled.bam"
    subsample_bam_for_qc "${bam}" "${temp_bam}" 5000000

    # Run RSeQC on subsampled BAM
    singularity exec ${RSEQC_CONTAINER} bam_stat.py -i ${temp_bam} > ${OUTPUT}/bamstat/${sample}.txt
    singularity exec ${RSEQC_CONTAINER} read_distribution.py -i ${temp_bam} -r ${BED} > ${OUTPUT}/read_distribution/${sample}.txt
    # ... etc

    # Clean up
    rm -f ${temp_bam} ${temp_bam}.bai
done
```

## Nextflow Implementation Considerations

For integration into the methyltna pipeline:

1. **Add subsampling process** before RSEQC_ANALYSIS subworkflow:
   ```groovy
   process SUBSAMPLE_BAM_FOR_QC {
       tag "$meta.id"

       input:
       tuple val(meta), path(bam), path(bai)
       val(max_reads)  // e.g., 5000000

       output:
       tuple val(meta), path("*_subsampled.bam"), path("*_subsampled.bam.bai")

       script:
       def total_reads = // count with samtools
       def fraction = max_reads / total_reads

       if (total_reads > max_reads) {
           """
           samtools view -s ${fraction} -b ${bam} > ${meta.id}_subsampled.bam
           samtools index ${meta.id}_subsampled.bam
           """
       } else {
           """
           ln -s ${bam} ${meta.id}_subsampled.bam
           ln -s ${bai} ${meta.id}_subsampled.bam.bai
           """
       }
   }
   ```

2. **Add parameter** to nextflow.config:
   ```groovy
   params {
       rseqc_max_reads = 5000000  // Maximum reads for RSeQC subsampling
   }
   ```

3. **Wire into workflow**:
   ```groovy
   // In workflows/methyltna.nf
   SUBSAMPLE_BAM_FOR_QC(
       SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai),
       params.rseqc_max_reads
   )

   RSEQC_ANALYSIS(
       SUBSAMPLE_BAM_FOR_QC.out.bam,
       SUBSAMPLE_BAM_FOR_QC.out.bai,
       PREPARE_REFERENCES.out.annotation_bed
   )
   ```

## Expected Performance Improvements

Based on real timing analysis from mot26 run (see `/home/marcus/runs/mot26/timing_analysis.txt`):

**Performance Characteristics by Sample Type:**
- **DNA samples (4M reads)**: ~24 min (~5.7 min per million reads)
- **TNA samples (35M reads avg)**: ~183 min (~5.4 min per million reads)
- **RNA samples (10M reads avg)**: ~86 min (~11.1 min per million reads)

**Key Insight**: RNA samples take **2x longer per million reads** than DNA/TNA due to better transcript coverage requiring more feature processing.

**Current (no subsampling):**
- 20 mixed samples = **~40 hours total**
- High variability: 24 min (DNA) to 204 min (TNA) per sample

**With 5M read subsampling:**
- DNA (4M): 24 min → no change (already <5M)
- TNA (35M): 183 min → 27 min (**6.8x speedup**)
- RNA (10M): 86 min → 56 min (**1.5x speedup**)
- **Total: ~10 hours for 20 samples (4x speedup!)**
- Much more manageable for production pipelines

**With 2M read subsampling:**
- DNA (4M): 24 min → 11 min (2.2x speedup)
- TNA (35M): 183 min → 11 min (**16.6x speedup!**)
- RNA (10M): 86 min → 22 min (3.9x speedup)
- **Total: ~4-5 hours for 20 samples (8-10x speedup!)**
- Recommended for production use

## Testing Plan

1. **Validate metrics consistency:**
   - Run RSeQC on full BAM vs. subsampled BAM for a few samples
   - Compare read_distribution percentages, gene body coverage curves
   - Should be nearly identical (within statistical noise)

2. **Test edge cases:**
   - Very low coverage samples (<1M reads)
   - Very high coverage samples (>50M reads)
   - Samples right at the threshold (5M reads)

3. **Performance benchmarking:**
   - Measure runtime with and without subsampling
   - Confirm expected speedup

4. **MultiQC integration:**
   - Ensure subsampled data displays correctly in MultiQC
   - Add note to report about subsampling strategy

## Additional Optimizations (Optional)

1. **Skip geneBody_coverage entirely** if runtime is still an issue
   - read_distribution already shows feature bias
   - Most valuable QC comes from read_distribution and inner_distance

2. **Reduce BED file complexity:**
   - Use filtered BED with only high-confidence protein-coding genes
   - E.g., 20K genes instead of 254K transcripts
   - Further ~10x speedup for geneBody_coverage

3. **Parallel processing:**
   - Already handled by Nextflow
   - Ensure adequate CPU allocation for samtools subsampling

## Implementation Priority

**Phase 1 (Recommended for immediate implementation):**
- Add subsampling to 5M reads before all RSeQC steps
- Test on a few samples
- Deploy to production

**Phase 2 (Future optimization):**
- Make threshold configurable via pipeline parameter
- Add subsampling stats to MultiQC report
- Consider per-tool subsampling (different thresholds for different tools)

## Questions for Implementer

1. Should subsampled BAMs be kept or deleted after RSeQC?
   - **Recommend**: Delete to save space (they're only for QC)

2. What should the default max_reads parameter be?
   - **Recommend**: 5M (good balance of speed and quality)
   - Make it configurable

3. Should this apply to all samples or only RNA samples?
   - **Recommend**: Only barcoded (RNA) samples
   - DNA samples are likely lower coverage anyway

4. Add to pipeline_info report?
   - **Recommend**: Yes, note which samples were subsampled and to what fraction

## Current Status

This document was created based on detailed timing analysis of the mot26 production run (see `/home/marcus/runs/mot26/timing_analysis.txt`), where RSeQC gene body coverage takes **~40 hours for 20 mixed samples** with extreme variability:
- DNA samples: ~24 minutes per sample
- TNA samples: ~183 minutes per sample (up to 3+ hours)
- RNA samples: ~86 minutes per sample

The 2x longer processing time per million reads for RNA samples (due to better transcript coverage) makes fixed-count subsampling particularly beneficial. The subsampling strategy has been validated in other RNA-seq pipelines (e.g., nf-core/rnaseq uses similar approaches) and is **essential for production viability**.

## References

- RSeQC documentation: http://rseqc.sourceforge.net/
- samtools subsampling: http://www.htslib.org/doc/samtools-view.html
- nf-core/rnaseq subsampling approach: https://github.com/nf-core/rnaseq

---

**Contact**: Pass this document to implementation team for integration into the methyltna pipeline.
