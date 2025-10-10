include { RNABARCODEEXTRACTION } from '../../../modules/local/rnabarcodeextraction/main'
include { RNA_STATS_SUMMARY     } from '../../../modules/local/rna_stats_summary/main'

workflow TNA_DECONVOLUTION {

    take:
    ch_reads_with_cfg  // channel: [ val(meta), [ fastq_1, fastq_2 ], path(config_file) ]

    main:

    ch_versions = Channel.empty()

    //
    // Extract RNA barcodes from paired-end FASTQ reads
    //
    RNABARCODEEXTRACTION ( ch_reads_with_cfg )
    ch_versions = ch_versions.mix(RNABARCODEEXTRACTION.out.versions.first())

    //
    // Create summary statistics from JSON reports for MultiQC
    //
    RNA_STATS_SUMMARY ( RNABARCODEEXTRACTION.out.json_report.map{ meta, json -> json }.collect() )
    ch_versions = ch_versions.mix(RNA_STATS_SUMMARY.out.versions)

    emit:
    barcoded_reads   = RNABARCODEEXTRACTION.out.barcoded_reads.map { meta, reads ->
        [meta + [id: "${meta.id}_barcoded", single_end: false], reads]
    }   // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    unbarcoded_reads = RNABARCODEEXTRACTION.out.unbarcoded_reads.map { meta, reads ->
        [meta + [id: "${meta.id}_unbarcoded", single_end: false], reads]
    } // channel: [ val(meta), [ fastq_1, fastq_2 ] ]
    json_report      = RNABARCODEEXTRACTION.out.json_report      // channel: [ val(meta), [ json ] ]
    txt_report       = RNABARCODEEXTRACTION.out.txt_report       // channel: [ val(meta), [ txt ] ]
    rna_stats_mqc    = RNA_STATS_SUMMARY.out.mqc                 // channel: [ rna_barcode_stats_mqc.txt ]

    versions         = ch_versions                               // channel: [ versions.yml ]
}
