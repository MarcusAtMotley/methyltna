//
// Read trimming workflow: FastQC -> TrimGalore -> FastQC
//

include { FASTQC                    } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_POST_TRIM } from '../../../modules/nf-core/fastqc/main'
include { TRIMGALORE } from '../../../modules/nf-core/trimgalore/main'

workflow READ_TRIMMING {

    take:
    ch_reads // channel: [ val(meta), [ fastq_1, fastq_2 ] ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: FastQC (pre-trim)
    //
    FASTQC (
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: TrimGalore
    //
    TRIMGALORE (
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())

    //
    // MODULE: FastQC (post-trim) - modify meta.id to avoid filename collision
    //
    ch_trimmed_reads_renamed = TRIMGALORE.out.reads.map { meta, reads ->
        def new_meta = meta.clone()
        new_meta.id = "${meta.id}_trimmed"
        [new_meta, reads]
    }
    
    FASTQC_POST_TRIM (
        ch_trimmed_reads_renamed
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST_TRIM.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_POST_TRIM.out.versions.first())

    emit:
    reads        = TRIMGALORE.out.reads       // channel: [ val(meta), [ fastq ] ]
    trimmed_reads = TRIMGALORE.out.reads      // channel: [ val(meta), [ fastq ] ] (alias for clarity)
    
    // QC outputs for MultiQC
    fastqc_raw_zip    = FASTQC.out.zip        // channel: [ val(meta), [ zip ] ]
    fastqc_raw_html   = FASTQC.out.html       // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip   = FASTQC_POST_TRIM.out.zip   // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html  = FASTQC_POST_TRIM.out.html  // channel: [ val(meta), [ html ] ]
    trim_log          = TRIMGALORE.out.log     // channel: [ val(meta), [ txt ] ]
    trim_html         = TRIMGALORE.out.html    // channel: [ val(meta), [ html ] ]
    trim_zip          = TRIMGALORE.out.zip     // channel: [ val(meta), [ zip ] ]
    
    multiqc_files     = ch_multiqc_files       // channel: [ files ]
    versions          = ch_versions            // channel: [ versions.yml ]
}