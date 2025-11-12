/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BCL_DEMULTIPLEX        } from '../subworkflows/nf-core/bcl_demultiplex/main'
include { READ_TRIMMING          } from '../subworkflows/local/read_trimming/main'
include { TNA_DECONVOLUTION      } from '../subworkflows/local/tna_deconvolution/main'
include { PREPARE_REFERENCES     } from '../subworkflows/local/prepare_references/main'
include { RSEQC_ANALYSIS         } from '../subworkflows/local/rseqc_analysis/main'
include { SUBSAMPLE_BAM_FOR_QC   } from '../modules/local/subsample_bam_for_qc/main'
include { STAR_ALIGN             } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_VIEW          } from '../modules/nf-core/samtools/view/main'
include { SAMTOOLS_SORT          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_FLAGSTAT      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS      } from '../modules/nf-core/samtools/idxstats/main'
include { SUBREAD_FEATURECOUNTS  } from '../modules/nf-core/subread/featurecounts/main'
include { BISCUIT_ALIGN          } from '../modules/nf-core/biscuit/align/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include { BISCUIT_PILEUP         } from '../modules/nf-core/biscuit/pileup/main'
include { BISCUIT_QC             } from '../modules/nf-core/biscuit/qc/main'
include { LOFREQ_VARIANT_CALLING as LOFREQ_RNA } from '../subworkflows/local/lofreq_variant_calling/main'
include { LOFREQ_SUMMARY         } from '../modules/local/lofreq_summary/main'
include { BISCUIT_METHYLATION_SUMMARY } from '../modules/local/biscuit_methylation_summary/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_methyltna_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METHYLTNA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Conditional logic: BCL demultiplexing or direct FASTQ input
    //
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

        // Optional sample filtering after demultiplexing
        ch_demux_fastq = BCL_DEMULTIPLEX.out.fastq

        // Apply regex filter if provided
        if (params.sample_filter_regex) {
            log.info "Filtering samples with regex: ${params.sample_filter_regex}"

            ch_demux_fastq = ch_demux_fastq
                .filter { meta, reads ->
                    meta.id =~ params.sample_filter_regex
                }
                .ifEmpty {
                    error """
                    ================================================================================
                    ERROR: Sample filter regex matched NO samples!

                    Regex pattern: ${params.sample_filter_regex}

                    Please check your --sample_filter_regex pattern.
                    Example patterns:
                      - ".*TNA.*"           # Matches samples containing "TNA"
                      - ".*[DTR]NA-EM.*"    # Matches DNA-EM, TNA-EM, RNA-EM
                      - "^CoB_.*"           # Matches samples starting with "CoB_"

                    Check demultiplexing reports to see actual sample names.
                    ================================================================================
                    """.stripIndent()
                }
        }

        // Apply sample list filter if provided
        if (params.sample_filter_file) {
            def sample_list = file(params.sample_filter_file, checkIfExists: true)
                .readLines()
                .collect { it.trim() }
                .findAll { it && !it.startsWith('#') }  // Remove empty lines and comments
            log.info "Filtering samples from file: ${params.sample_filter_file} (${sample_list.size()} samples in filter)"

            ch_demux_fastq = ch_demux_fastq.filter { meta, reads ->
                meta.id in sample_list
            }
            .ifEmpty {
                error """
                ================================================================================
                ERROR: Sample filter file matched NO samples!

                Filter file: ${params.sample_filter_file}
                Samples in filter: ${sample_list.join(', ')}

                None of the samples in your filter file were found in the demultiplexed data.
                Please check that sample IDs match exactly (case-sensitive).
                ================================================================================
                """.stripIndent()
            }
        }

        ch_fastq_for_trimming = ch_demux_fastq
        ch_versions = ch_versions.mix(BCL_DEMULTIPLEX.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BCL_DEMULTIPLEX.out.reports.collect{it[1]})

    } else {
        // Direct FASTQ input path
        ch_fastq_for_trimming = ch_samplesheet
    }

    //
    // SUBWORKFLOW: Read trimming (FastQC → TrimGalore → FastQC)
    //
    READ_TRIMMING(ch_fastq_for_trimming)
    ch_versions = ch_versions.mix(READ_TRIMMING.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(READ_TRIMMING.out.multiqc_files)

    //
    // SUBWORKFLOW: RNA barcode deconvolution
    //
    final cfg_file = file(params.rna_barcode_config, checkIfExists: true)

    ch_reads_with_cfg = READ_TRIMMING.out.reads.map { meta, reads ->
        tuple(meta, reads, cfg_file)
    }

    TNA_DECONVOLUTION(ch_reads_with_cfg)
    ch_versions = ch_versions.mix(TNA_DECONVOLUTION.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.txt_report.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.json_report.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(TNA_DECONVOLUTION.out.rna_stats_mqc)

    //
    // SUBWORKFLOW: Prepare references (download and index if needed)
    //
    PREPARE_REFERENCES()
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)

    // Convert reference channels to value channels for multiple consumption
    ch_star_index = PREPARE_REFERENCES.out.star_index.first()
    ch_gtf = PREPARE_REFERENCES.out.annotation_gtf.first()
    ch_biscuit_index = PREPARE_REFERENCES.out.biscuit_index.first()
    ch_genome_fasta = PREPARE_REFERENCES.out.genome_fasta.first()
    ch_genome_fai = PREPARE_REFERENCES.out.genome_fai.first()

    //
    // RNA PATH: Use barcoded reads (includes all RNA regardless of single/double tagging)
    //
    ch_rna_reads = TNA_DECONVOLUTION.out.barcoded_reads

    //
    // MODULE: STAR alignment for RNA reads (splice-aware)
    //
    STAR_ALIGN(
        ch_rna_reads,
        ch_star_index,
        ch_gtf,
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
        STAR_ALIGN.out.sam.map { meta, sam -> [meta, sam, []] },
        [[],[]],
        [],
        ''
    )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions)

    //
    // MODULE: Sort BAM files for downstream analysis
    //
    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam,
        [[],[]],  // No reference fasta needed for coordinate sorting
        'bai'     // Generate BAI index format
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    //
    // MODULE: Optionally subsample BAMs for faster RSeQC QC
    //
    if (!params.skip_rseqc_subsampling && params.rseqc_subsample_reads > 0) {
        // Subsample BAMs to improve RSeQC performance
        ch_rseqc_bam_input = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_SORT.out.bai)

        SUBSAMPLE_BAM_FOR_QC(
            ch_rseqc_bam_input,
            params.rseqc_subsample_reads
        )
        ch_versions = ch_versions.mix(SUBSAMPLE_BAM_FOR_QC.out.versions)

        // Use subsampled BAMs for RSeQC
        ch_rseqc_bam = SUBSAMPLE_BAM_FOR_QC.out.bam.map{ meta, bam, bai -> [meta, bam] }
        ch_rseqc_bai = SUBSAMPLE_BAM_FOR_QC.out.bam.map{ meta, bam, bai -> [meta, bai] }
    } else {
        // Use full-depth BAMs for RSeQC (no subsampling)
        ch_rseqc_bam = SAMTOOLS_SORT.out.bam
        ch_rseqc_bai = SAMTOOLS_SORT.out.bai
    }

    //
    // SUBWORKFLOW: RSeQC RNA-seq Quality Control
    //
    RSEQC_ANALYSIS(
        ch_rseqc_bam,
        ch_rseqc_bai,
        PREPARE_REFERENCES.out.annotation_bed
    )
    ch_versions = ch_versions.mix(RSEQC_ANALYSIS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_ANALYSIS.out.multiqc_files)

    //
    // MODULE: FeatureCounts for gene quantification
    //
    ch_gtf_file = PREPARE_REFERENCES.out.annotation_gtf.map{ it[1] }.first()

    SUBREAD_FEATURECOUNTS(
        SAMTOOLS_SORT.out.bam.combine(ch_gtf_file)
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]})

    //
    // SUBWORKFLOW: LoFreq variant calling on RNA reads
    //
    LOFREQ_RNA(
        SAMTOOLS_SORT.out.bam,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(LOFREQ_RNA.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.flagstat.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_RNA.out.idxstats.collect{it[1]})

    //
    // MODULE: Generate LoFreq variant summary for MultiQC
    //
    LOFREQ_SUMMARY(
        LOFREQ_RNA.out.vcf.collect{it[1]}
    )
    ch_versions = ch_versions.mix(LOFREQ_SUMMARY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(LOFREQ_SUMMARY.out.json)

    //
    // DNA PATH: Unbarcoded reads
    //
    ch_dna_reads = TNA_DECONVOLUTION.out.unbarcoded_reads

    //
    // MODULE: Biscuit alignment for DNA reads (EM-seq)
    //
    BISCUIT_ALIGN(
        ch_dna_reads,
        ch_genome_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_ALIGN.out.versions)

    //
    // MODULE: Mark duplicates in aligned DNA reads and create BAI index
    //
    // Note: Picard generates BAI automatically with --CREATE_INDEX flag
    PICARD_MARKDUPLICATES(
        BISCUIT_ALIGN.out.bam,
        [[],[]],  // No reference fasta needed
        [[],[]]   // No reference fasta index needed
    )
    ch_versions = ch_versions.mix(PICARD_MARKDUPLICATES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{it[1]})

    //
    // MODULE: SAMtools flagstat for DNA alignment statistics
    //
    SAMTOOLS_FLAGSTAT(
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, by: [0])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})

    //
    // MODULE: SAMtools idxstats for DNA alignment statistics
    //
    SAMTOOLS_IDXSTATS(
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, by: [0])
    )
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]})

    //
    // MODULE: Mosdepth coverage analysis for DNA reads
    //
    MOSDEPTH(
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, by: [0]).map { meta, bam, bai -> [meta, bam, bai, []] },
        ch_genome_fasta
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})

    //
    // MODULE: Biscuit pileup for methylation calling and variant detection
    //
    // Note: BISCUIT_PILEUP expects (normal_bams, normal_bais, tumor_bam, tumor_bai)
    // For germline: provide BAM as normal_bams and leave tumor channels empty
    BISCUIT_PILEUP(
        PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, by: [0]).map { meta, bam, bai -> [meta, [bam], [bai], [], []] },
        ch_genome_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_PILEUP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BISCUIT_PILEUP.out.vcf.collect{it[1]})

    //
    // MODULE: Biscuit QC for methylation quality metrics
    //
    // Note: EM-seq data may show -nan for conversion rates since EM-seq uses enzymatic
    // conversion (not bisulfite), so traditional conversion efficiency doesn't apply.
    // Actual methylation data is in the BISCUIT_PILEUP VCF files.
    BISCUIT_QC(
        PICARD_MARKDUPLICATES.out.bam,
        ch_genome_fasta,
        ch_biscuit_index
    )
    ch_versions = ch_versions.mix(BISCUIT_QC.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BISCUIT_QC.out.reports.collect{it[1]})

    //
    // MODULE: Extract methylation summary statistics from Biscuit VCF files
    //
    BISCUIT_METHYLATION_SUMMARY(
        BISCUIT_PILEUP.out.vcf.collect{it[1]}
    )
    ch_versions = ch_versions.mix(BISCUIT_METHYLATION_SUMMARY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(BISCUIT_METHYLATION_SUMMARY.out.json)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'methyltna_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
