/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { BCL_DEMULTIPLEX        } from '../subworkflows/nf-core/bcl_demultiplex/main'
include { METHYLSNP_HAIRPIN_PROCESSING } from '../subworkflows/local/methylsnp_hairpin_processing/main'
include { HAIRPIN_RESOLUTION_STATS_SUMMARY } from '../modules/local/hairpin_resolution_stats_summary/main'
include { TNA_DECONVOLUTION      } from '../subworkflows/local/tna_deconvolution/main'
include { PREPARE_REFERENCES     } from '../subworkflows/local/prepare_references/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { METHYLSNP_PROCESSING   } from '../modules/local/methylsnp_processing/main'
include { METHYLSNP_EXTRACTION   } from '../modules/local/methylsnp_extraction/main'
include { LOFREQ_VARIANT_CALLING as LOFREQ_RNA } from '../subworkflows/local/lofreq_variant_calling/main'
include { LOFREQ_VARIANT_CALLING as LOFREQ_DNA } from '../subworkflows/local/lofreq_variant_calling/main'
include { LOFREQ_VARIANT_CALLING as LOFREQ_SINGLE } from '../subworkflows/local/lofreq_variant_calling/main'
include { LOFREQ_SUMMARY         } from '../modules/local/lofreq_summary/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { SUBREAD_FEATURECOUNTS  } from '../modules/nf-core/subread/featurecounts/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_modulestesting_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MODULESTESTING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_fastq_for_qc = Channel.empty()

    // Conditional logic: BCL demultiplexing or direct FASTQ input
    if (params.bcl_input_dir && params.bcl_samplesheet) {
        // BCL demultiplexing path
        
        // Create flowcell channel for BCL_DEMULTIPLEX
        // Format: [[id:"", lane:""], samplesheet.csv, path/to/bcl/files]
        ch_flowcell = Channel.of([
            [id: "flowcell", lane: ""],
            file(params.bcl_samplesheet, checkIfExists: true),
            file(params.bcl_input_dir, checkIfExists: true)
        ])

        BCL_DEMULTIPLEX(
            ch_flowcell,
            params.demultiplexer
        )
        
        ch_fastq_for_qc = BCL_DEMULTIPLEX.out.fastq
        ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCL_DEMULTIPLEX.out.reports.collect{it[1]})
        
    } else {
        // Direct FASTQ input path
        ch_fastq_for_qc = ch_samplesheet
    }

    //
    // MODULE: Optional FastQC on raw reads (pre-processing QC only)
    //
    if (!params.skip_fastqc) {
        FASTQC(ch_fastq_for_qc)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
        ch_versions = ch_versions.mix(FASTQC.out.versions)
    }

    //
    // SUBWORKFLOW: Hairpin processing FIRST (v2 protocol - works on all reads)
    //
    if (!params.skip_methylsnp_analysis) {
        METHYLSNP_HAIRPIN_PROCESSING(ch_fastq_for_qc)
        ch_versions = ch_versions.mix(METHYLSNP_HAIRPIN_PROCESSING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_HAIRPIN_PROCESSING.out.illumina_reports.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_HAIRPIN_PROCESSING.out.hairpin_reports.collect{it[1]})

        ch_processed_reads = METHYLSNP_HAIRPIN_PROCESSING.out.resolved_reads
        ch_methylation_reports = METHYLSNP_HAIRPIN_PROCESSING.out.methylation_report

        // Collect and summarize hairpin resolution statistics
        HAIRPIN_RESOLUTION_STATS_SUMMARY(
            METHYLSNP_HAIRPIN_PROCESSING.out.hairpin_resolution_stats.collect{it[1]}
        )
        ch_versions = ch_versions.mix(HAIRPIN_RESOLUTION_STATS_SUMMARY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(HAIRPIN_RESOLUTION_STATS_SUMMARY.out.mqc_yaml)
    } else {
        // Skip hairpin processing - use raw reads
        ch_processed_reads = ch_fastq_for_qc
    }

    //
    // SUBWORKFLOW: RNA barcode deconvolution (AFTER hairpin processing)
    //
    if (!params.skip_rna_deconvolution && params.rna_barcode_config) {
        // Make sure the file gets staged
        final cfg_file = file(params.rna_barcode_config, checkIfExists: true)

        // Attach cfg to every sample tuple: (meta, reads, cfg)
        def ch_reads_with_cfg = ch_processed_reads.map { meta, reads ->
            tuple(meta, reads, cfg_file)
        }

        // Call the subworkflow with a single channel of 3-tuples
        TNA_DECONVOLUTION(ch_reads_with_cfg)

        ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.txt_report.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.json_report.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.rna_stats_mqc)
        ch_versions = ch_versions.mix(TNA_DECONVOLUTION.out.versions)

        ch_rna_reads = TNA_DECONVOLUTION.out.barcoded_reads
        ch_dna_reads = TNA_DECONVOLUTION.out.unbarcoded_reads

        // Duplicate methylation reports with modified IDs to match barcoded/unbarcoded reads
        // Must use flatMap to create both variants from each input item (channels can only be consumed once)
        // IMPORTANT: Must match ALL meta fields from TNA_DECONVOLUTION (id suffix AND single_end: true)
        ch_methylation_reports_split = ch_methylation_reports
            .flatMap { meta, report ->
                [
                    [meta + [id: "${meta.id}_barcoded", single_end: true], report],
                    [meta + [id: "${meta.id}_unbarcoded", single_end: true], report]
                ]
            }
    } else {
        // No RNA barcode extraction - all reads treated as DNA
        ch_dna_reads = ch_processed_reads
        ch_methylation_reports_split = ch_methylation_reports
    }

    //
    // SUBWORKFLOW: Prepare references (download and index if needed)
    //
    if (!params.skip_methylsnp_analysis) {
        PREPARE_REFERENCES()
        ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

        //
        // CONDITIONAL ALIGNMENT: Dual or single processing
        //
        if (!params.skip_rna_deconvolution && params.rna_barcode_config) {
            // Dual processing: RNA → STAR, DNA → Bowtie2

            //
            // MODULE: STAR alignment for RNA reads (splice-aware)
            //
            // Convert reference channels to value channels for multiple consumption
            ch_star_index_val = PREPARE_REFERENCES.out.star_index.first()
            ch_gtf_val = PREPARE_REFERENCES.out.annotation_gtf.first()
            ch_bowtie2_index_val = PREPARE_REFERENCES.out.bowtie2_index.first()
            ch_genome_fasta_val = PREPARE_REFERENCES.out.genome_fasta.first()

            STAR_ALIGN(
                ch_rna_reads,
                ch_star_index_val,
                ch_gtf_val,
                false,          // star_ignore_sjdbgtf
                'ILLUMINA',     // seq_platform
                ''              // seq_center
            )
            ch_versions = ch_versions.mix(STAR_ALIGN.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(STAR_ALIGN.out.log_final.collect{it[1]})

            //
            // MODULE: Convert SAM to BAM for FeatureCounts
            //
            SAMTOOLS_VIEW(
                STAR_ALIGN.out.sam.map { meta, sam -> [meta, sam, []] },  // Add empty index
                [[],[]],                                                     // No reference needed
                [],                                                          // No qname list
                ''                                                           // No index format
            )
            ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

            //
            // MODULE: FeatureCounts for gene quantification
            //
            // Extract GTF file only (not meta map) for featureCounts
            ch_gtf_file = PREPARE_REFERENCES.out.annotation_gtf.map{ it[1] }.first()

            SUBREAD_FEATURECOUNTS(
                SAMTOOLS_VIEW.out.bam.combine(ch_gtf_file)
            )
            ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})

            //
            // MODULE: Bowtie2 alignment for DNA reads (genomic)
            //
            BOWTIE2_ALIGN(
                ch_dna_reads,
                ch_bowtie2_index_val,
                ch_genome_fasta_val,
                false,          // save_unaligned
                false           // sort_bam
            )
            ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})

            //
            // Combine RNA and DNA alignments for processing
            //
            ch_all_sam = STAR_ALIGN.out.sam.mix(BOWTIE2_ALIGN.out.sam)

            // Combine SAM with methylation reports (using split reports with matching IDs)
            ch_combined_for_processing = ch_all_sam
                .join(ch_methylation_reports_split, by: [0])
                .combine(ch_genome_fasta_val)

            //
            // MODULE: SAM processing (mark unique, duplicates, add XM tags)
            //
            METHYLSNP_PROCESSING(
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta, sam]
                },
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta, report]
                },
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta3, fasta]
                }
            )
            ch_versions = ch_versions.mix(METHYLSNP_PROCESSING.out.versions)

            //
            // MODULE: Bismark methylation extraction
            //
            METHYLSNP_EXTRACTION(METHYLSNP_PROCESSING.out.sorted_bam)
            ch_versions = ch_versions.mix(METHYLSNP_EXTRACTION.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_EXTRACTION.out.report.collect{it[1]})

            // Split processed BAMs back into RNA and DNA for variant calling
            ch_rna_bam = METHYLSNP_PROCESSING.out.sorted_bam.filter { meta, bam -> meta.id.contains('_barcoded') }
            ch_dna_bam = METHYLSNP_PROCESSING.out.sorted_bam.filter { meta, bam -> meta.id.contains('_unbarcoded') }

            //
            // Variant calling on RNA reads
            //
            ch_genome_fai_val = PREPARE_REFERENCES.out.genome_fai.first()

            LOFREQ_RNA (
                ch_rna_bam,
                ch_genome_fasta_val,
                ch_genome_fai_val
            )
            ch_versions = ch_versions.mix(LOFREQ_RNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.vcf.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.flagstat.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.idxstats.collect{it[1]})

            //
            // Variant calling on DNA reads
            //
            LOFREQ_DNA (
                ch_dna_bam,
                ch_genome_fasta_val,
                ch_genome_fai_val
            )
            ch_versions = ch_versions.mix(LOFREQ_DNA.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_DNA.out.vcf.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_DNA.out.flagstat.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_DNA.out.idxstats.collect{it[1]})

            // Collect VCF files from RNA and DNA for summary
            ch_all_vcfs = LOFREQ_RNA.out.vcf.collect{it[1]}.mix(LOFREQ_DNA.out.vcf.collect{it[1]}).collect()

        } else {
            // Single processing: all reads treated as DNA → Bowtie2

            // Convert reference channels to value channels for multiple consumption
            ch_bowtie2_index_val = PREPARE_REFERENCES.out.bowtie2_index.first()
            ch_genome_fasta_val = PREPARE_REFERENCES.out.genome_fasta.first()
            ch_genome_fai_val = PREPARE_REFERENCES.out.genome_fai.first()

            //
            // MODULE: Bowtie2 alignment for all reads (genomic)
            //
            BOWTIE2_ALIGN(
                ch_dna_reads,
                ch_bowtie2_index_val,
                ch_genome_fasta_val,
                false,          // save_unaligned
                false           // sort_bam
            )
            ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log.collect{it[1]})

            // Combine SAM with methylation reports (using split reports)
            ch_combined_for_processing = BOWTIE2_ALIGN.out.sam
                .join(ch_methylation_reports_split, by: [0])
                .combine(ch_genome_fasta_val)

            //
            // MODULE: SAM processing (mark unique, duplicates, add XM tags)
            //
            METHYLSNP_PROCESSING(
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta, sam]
                },
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta, report]
                },
                ch_combined_for_processing.map { meta, sam, report, meta3, fasta ->
                    [meta3, fasta]
                }
            )
            ch_versions = ch_versions.mix(METHYLSNP_PROCESSING.out.versions)

            //
            // MODULE: Bismark methylation extraction
            //
            METHYLSNP_EXTRACTION(METHYLSNP_PROCESSING.out.sorted_bam)
            ch_versions = ch_versions.mix(METHYLSNP_EXTRACTION.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(METHYLSNP_EXTRACTION.out.report.collect{it[1]})

            //
            // Variant calling on all reads
            //
            LOFREQ_SINGLE (
                METHYLSNP_PROCESSING.out.sorted_bam,
                ch_genome_fasta_val,
                ch_genome_fai_val
            )
            ch_versions = ch_versions.mix(LOFREQ_SINGLE.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SINGLE.out.vcf.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SINGLE.out.flagstat.collect{it[1]})
            ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SINGLE.out.idxstats.collect{it[1]})

            // Collect all VCF files for summary
            ch_all_vcfs = LOFREQ_SINGLE.out.vcf.collect{it[1]}
        }

        // Generate LoFreq summary statistics for MultiQC
        LOFREQ_SUMMARY(ch_all_vcfs)
        ch_versions = ch_versions.mix(LOFREQ_SUMMARY.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SUMMARY.out.json)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'modulestesting_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
