#!/usr/bin/env nextflow

//
// LOFREQ_VARIANT_CALLING: Low-frequency variant calling using LoFreq
//

include { LOFREQ_CALLPARALLEL } from '../../../modules/nf-core/lofreq/callparallel/main'
include { SAMTOOLS_INDEX      } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS   } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT   } from '../../../modules/nf-core/samtools/flagstat/main'

workflow LOFREQ_VARIANT_CALLING {

    take:
    ch_bam_files       // channel: [meta, bam]
    ch_reference_fasta // channel: [meta, fasta]
    ch_reference_fai   // channel: [meta, fai]

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Index BAM files for variant calling
    //
    SAMTOOLS_INDEX(ch_bam_files)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Check BAM files for reads using idxstats
    ch_bam_with_index = ch_bam_files
        .join(SAMTOOLS_INDEX.out.bai, by: [0])

    SAMTOOLS_IDXSTATS(ch_bam_with_index)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    // Generate flagstat reports for MultiQC
    SAMTOOLS_FLAGSTAT(ch_bam_with_index)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

    // Filter out empty BAM files by checking idxstats output
    ch_non_empty_bams = ch_bam_with_index
        .join(SAMTOOLS_IDXSTATS.out.idxstats, by: [0])
        .filter { meta, bam, bai, idxstats ->
            // Read the idxstats file to check for mapped reads
            def statsText = idxstats.text.trim()
            def totalReads = 0

            statsText.split('\n').each { line ->
                if (line && !line.startsWith('*')) {
                    def parts = line.split('\t')
                    if (parts.size() >= 3) {
                        totalReads += parts[2] as Integer
                    }
                }
            }

            if (totalReads == 0) {
                log.warn "Skipping LoFreq variant calling for ${meta.id}: BAM file contains no mapped reads"
                return false
            }
            return true
        }
        .map { meta, bam, bai, idxstats -> [meta, bam, bai] }

    // Prepare channels for LoFreq - need to combine with reference files
    ch_lofreq_input = ch_non_empty_bams
        .map { meta, bam, bai -> [meta, bam, bai, []] }  // Add empty intervals
        .combine(ch_reference_fasta)
        .combine(ch_reference_fai)

    //
    // MODULE: LoFreq variant calling
    //
    LOFREQ_CALLPARALLEL(
        ch_lofreq_input.map { meta, bam, bai, intervals, meta2, fasta, meta3, fai ->
            [meta, bam, bai, intervals]
        },
        ch_lofreq_input.map { meta, bam, bai, intervals, meta2, fasta, meta3, fai ->
            [meta2, fasta]
        },
        ch_lofreq_input.map { meta, bam, bai, intervals, meta2, fasta, meta3, fai ->
            [meta3, fai]
        }
    )
    ch_versions = ch_versions.mix(LOFREQ_CALLPARALLEL.out.versions)

    emit:
    vcf      = LOFREQ_CALLPARALLEL.out.vcf    // channel: [meta, vcf.gz]
    tbi      = LOFREQ_CALLPARALLEL.out.tbi    // channel: [meta, vcf.gz.tbi]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [meta, flagstat]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [meta, idxstats]
    versions = ch_versions                    // channel: [versions.yml]
}