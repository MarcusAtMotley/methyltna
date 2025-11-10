/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: RSeQC RNA-seq Quality Control Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Comprehensive RNA-seq QC using RSeQC toolkit
    - BAM statistics
    - Read distribution across genomic features
    - Inner distance analysis for paired-end data
    - Gene body coverage analysis (custom module)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RSEQC_BAMSTAT            } from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_READDISTRIBUTION   } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_INNERDISTANCE      } from '../../../modules/nf-core/rseqc/innerdistance/main'
include { RSEQC_GENEBODYCOVERAGE   } from '../../../modules/local/rseqc/genebodycoverage/main'

workflow RSEQC_ANALYSIS {
    take:
    bam         // channel: [ val(meta), path(bam) ]
    bed         // channel: [ val(meta), path(bed) ] - Pre-converted BED file from PREPARE_REFERENCES

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Extract BED file from channel (remove meta)
    ch_bed = bed.map { meta, bed_file -> bed_file }

    //
    // MODULE: BAM statistics
    //
    RSEQC_BAMSTAT(
        bam
    )
    ch_versions = ch_versions.mix(RSEQC_BAMSTAT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.collect{it[1]})

    //
    // MODULE: Read distribution across genomic features
    //
    RSEQC_READDISTRIBUTION(
        bam,
        ch_bed
    )
    ch_versions = ch_versions.mix(RSEQC_READDISTRIBUTION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_READDISTRIBUTION.out.txt.collect{it[1]})

    //
    // MODULE: Inner distance analysis (paired-end only)
    //
    RSEQC_INNERDISTANCE(
        bam,
        ch_bed
    )
    ch_versions = ch_versions.mix(RSEQC_INNERDISTANCE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_INNERDISTANCE.out.distance.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_INNERDISTANCE.out.freq.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_INNERDISTANCE.out.mean.collect{it[1]})

    //
    // MODULE: Gene body coverage analysis (custom)
    //
    RSEQC_GENEBODYCOVERAGE(
        bam,
        ch_bed
    )
    ch_versions = ch_versions.mix(RSEQC_GENEBODYCOVERAGE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_GENEBODYCOVERAGE.out.txt.collect{it[1]})

    emit:
    bamstat_txt          = RSEQC_BAMSTAT.out.txt              // channel: [ val(meta), path(txt) ]
    readdist_txt         = RSEQC_READDISTRIBUTION.out.txt     // channel: [ val(meta), path(txt) ]
    innerdist_distance   = RSEQC_INNERDISTANCE.out.distance   // channel: [ val(meta), path(txt) ]
    innerdist_freq       = RSEQC_INNERDISTANCE.out.freq       // channel: [ val(meta), path(txt) ]
    innerdist_mean       = RSEQC_INNERDISTANCE.out.mean       // channel: [ val(meta), path(txt) ]
    innerdist_pdf        = RSEQC_INNERDISTANCE.out.pdf        // channel: [ val(meta), path(pdf) ]
    genebody_txt         = RSEQC_GENEBODYCOVERAGE.out.txt     // channel: [ val(meta), path(txt) ]
    genebody_pdf         = RSEQC_GENEBODYCOVERAGE.out.pdf     // channel: [ val(meta), path(pdf) ]

    multiqc_files        = ch_multiqc_files                   // channel: path(*)
    versions             = ch_versions                        // channel: path(versions.yml)
}
