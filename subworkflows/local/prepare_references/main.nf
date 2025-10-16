#!/usr/bin/env nextflow

//
// PREPARE_REFERENCES: Smart reference download and indexing with caching
//

include { DOWNLOAD_REFERENCES } from '../../../modules/local/download_references/main'
include { DOWNLOAD_CLOUD_CACHE } from '../../../modules/local/download_cloud_cache/main'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'
include { BOWTIE2_BUILD } from '../../../modules/nf-core/bowtie2/build/main'
include { BISCUIT_INDEX } from '../../../modules/nf-core/biscuit/index/main'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME } from '../../../modules/nf-core/samtools/faidx/main'

workflow PREPARE_REFERENCES {

    take:
    // No inputs - references determined by parameters

    main:

    ch_versions = Channel.empty()

    // Extract filenames from URLs for cache checking (genome-agnostic)
    def genome_filename = params.genome_fasta.tokenize('/').last()
    def gtf_filename = params.annotation_gtf.tokenize('/').last()
    def genome_id = genome_filename.replaceAll(/\.fa(sta)?$/, '')

    // Check for locally cached files
    def genome_cached = file("${params.reference_cache_dir}/fasta/${genome_filename}").exists()
    def gtf_cached = file("${params.reference_cache_dir}/${gtf_filename}").exists()

    // Determine if we need to download based on cache status, provided paths, and force flags
    def need_download = false
    def need_cloud_download = false

    if (params.genome_fasta.startsWith('gs://') && (!genome_cached || params.force_redownload_references)) {
        // Check if available in cloud cache before downloading from original source
        if (!params.force_redownload_references && params.cloud_reference_cache) {
            need_cloud_download = true
        } else {
            need_download = true
        }
    }
    if (params.annotation_gtf.startsWith('gs://') && (!gtf_cached || params.force_redownload_references)) {
        // Check if available in cloud cache before downloading from original source
        if (!params.force_redownload_references && params.cloud_reference_cache) {
            need_cloud_download = true
        } else {
            need_download = true
        }
    }

    if (need_cloud_download) {
        //
        // MODULE: Check cloud cache and download if available (FASTA, GTF, and all indexes)
        //
        DOWNLOAD_CLOUD_CACHE(genome_filename, gtf_filename, genome_id)
        ch_versions = ch_versions.mix(DOWNLOAD_CLOUD_CACHE.out.versions)

        // Use cloud cache outputs if available, otherwise will fall back to original source
        ch_genome_fasta = DOWNLOAD_CLOUD_CACHE.out.genome_fasta
            .ifEmpty{ file("${params.reference_cache_dir}/fasta/${genome_filename}") }
            .map { fasta -> [[id:'genome'], fasta] }

        ch_annotation_gtf = DOWNLOAD_CLOUD_CACHE.out.gtf
            .ifEmpty{ file("${params.reference_cache_dir}/${gtf_filename}") }
            .map { gtf -> [[id:'annotation'], gtf] }

        // Check if cloud download succeeded
        genome_cached = file("${params.reference_cache_dir}/fasta/${genome_filename}").exists()
        gtf_cached = file("${params.reference_cache_dir}/${gtf_filename}").exists()

        // If still not cached after cloud attempt, need to download from original source
        if (!genome_cached || !gtf_cached) {
            need_download = true
        }
    }

    if (need_download) {
        //
        // MODULE: Download references from original GCS source if not in cache
        //
        DOWNLOAD_REFERENCES()
        ch_versions = ch_versions.mix(DOWNLOAD_REFERENCES.out.versions)

        ch_genome_fasta = DOWNLOAD_REFERENCES.out.genome_fasta.map { fasta -> [[id:'genome'], fasta] }
        ch_annotation_gtf = DOWNLOAD_REFERENCES.out.gtf.map { gtf -> [[id:'annotation'], gtf] }
    } else {
        // Use cached files if available, otherwise use provided local paths
        if (genome_cached && params.genome_fasta.startsWith('gs://')) {
            ch_genome_fasta = Channel.fromPath("${params.reference_cache_dir}/fasta/${genome_filename}")
                .map { fasta -> [[id:'genome'], fasta] }
        } else {
            ch_genome_fasta = Channel.fromPath(params.genome_fasta)
                .map { fasta -> [[id:'genome'], fasta] }
        }

        if (gtf_cached && params.annotation_gtf.startsWith('gs://')) {
            ch_annotation_gtf = Channel.fromPath("${params.reference_cache_dir}/${gtf_filename}")
                .map { gtf -> [[id:'annotation'], gtf] }
        } else {
            ch_annotation_gtf = Channel.fromPath(params.annotation_gtf)
                .map { gtf -> [[id:'annotation'], gtf] }
        }
    }

    // Build or use STAR genome index - check local cache (which may have been populated by cloud download)
    def star_index_dir = "${params.reference_cache_dir}/star_indexes/${genome_id}"
    def star_index_exists = file("${star_index_dir}").exists() &&
                            file("${star_index_dir}/Genome").exists()

    if (params.star_index != null) {
        // User provided explicit STAR index path
        ch_star_index = Channel.fromPath(params.star_index)
            .map { dir -> [[id:'genome'], dir] }
    } else if (star_index_exists && !params.force_rebuild_indexes) {
        // Use locally cached STAR index
        ch_star_index = Channel.fromPath(star_index_dir)
            .map { dir -> [[id:'genome'], dir] }
    } else {
        //
        // MODULE: Build STAR genome index (cloud cache not available or download failed)
        //
        // Add genome_id to meta for publishDir path
        def ch_genome_with_id = ch_genome_fasta.map { meta, fasta ->
            [meta + [genome_id: genome_id], fasta]
        }
        def ch_gtf_with_id = ch_annotation_gtf.map { meta, gtf ->
            [meta + [genome_id: genome_id], gtf]
        }
        STAR_GENOMEGENERATE(ch_genome_with_id, ch_gtf_with_id)
        ch_star_index = STAR_GENOMEGENERATE.out.index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }


    // Build or use Bowtie2 genome index - check local cache (which may have been populated by cloud download)
    def bowtie2_index_dir = "${params.reference_cache_dir}/bowtie2_indexes/${genome_id}"
    def bowtie2_index_exists = file("${bowtie2_index_dir}").exists() &&
                                file("${bowtie2_index_dir}").list().any { it.endsWith('.bt2') }

    if (params.bowtie2_index != null) {
        // User provided explicit Bowtie2 index path
        ch_bowtie2_index = Channel.fromPath(params.bowtie2_index)
            .map { dir -> [[id:'genome'], dir] }
    } else if (bowtie2_index_exists && !params.force_rebuild_indexes) {
        // Use cached Bowtie2 index
        ch_bowtie2_index = Channel.fromPath(bowtie2_index_dir)
            .map { dir -> [[id:'genome'], dir] }
    } else {
        //
        // MODULE: Build Bowtie2 genome index
        //
        BOWTIE2_BUILD(ch_genome_fasta)
        ch_bowtie2_index = BOWTIE2_BUILD.out.index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }


    // Build or use Biscuit genome index - check local cache (which may have been populated by cloud download)
    def biscuit_index_dir = "${params.reference_cache_dir}/biscuit_indexes/${genome_id}"
    def biscuit_index_exists = file("${biscuit_index_dir}").exists() &&
                                file("${biscuit_index_dir}/${genome_filename}.dau.bwt").exists()

    if (params.biscuit_index != null) {
        // User provided explicit Biscuit index path
        ch_biscuit_index = Channel.fromPath(params.biscuit_index)
            .map { dir -> [[id:'genome'], dir] }
    } else if (biscuit_index_exists && !params.force_rebuild_indexes) {
        // Use cached Biscuit index
        ch_biscuit_index = Channel.fromPath(biscuit_index_dir)
            .map { dir -> [[id:'genome'], dir] }
    } else {
        //
        // MODULE: Build Biscuit genome index
        //
        // Add genome_id to meta for publishDir path
        def ch_genome_with_id = ch_genome_fasta.map { meta, fasta ->
            [meta + [genome_id: genome_id], fasta]
        }
        BISCUIT_INDEX(ch_genome_with_id)
        ch_biscuit_index = BISCUIT_INDEX.out.index
        ch_versions = ch_versions.mix(BISCUIT_INDEX.out.versions)
    }


    //
    // MODULE: Index FASTA files for variant calling (with caching)
    //
    // Check for cached FAI file
    def fai_file_path = genome_cached ? "${params.reference_cache_dir}/fasta/${genome_filename}.fai" : null
    def fai_exists = fai_file_path ? file(fai_file_path).exists() : false

    if (fai_exists && !params.force_rebuild_indexes) {
        // Use cached FAI file
        ch_genome_fai = ch_genome_fasta.combine(Channel.fromPath(fai_file_path))
            .map { meta_fasta, fasta, fai -> [meta_fasta, fai] }
    } else {
        // Generate FAI file
        SAMTOOLS_FAIDX_GENOME(
            ch_genome_fasta,
            [[],[]],  // No existing FAI file
            false     // Don't generate sizes file
        )
        ch_genome_fai = SAMTOOLS_FAIDX_GENOME.out.fai
        ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions)
    }


    emit:
    genome_fasta        = ch_genome_fasta              // channel: [meta, genome.fa]
    annotation_gtf      = ch_annotation_gtf            // channel: [meta, annotation.gtf]
    star_index          = ch_star_index                // channel: [meta, star_index_dir]
    bowtie2_index       = ch_bowtie2_index             // channel: [meta, bowtie2_index_dir]
    biscuit_index       = ch_biscuit_index             // channel: [meta, biscuit_index_dir]
    genome_fai          = ch_genome_fai                // channel: [meta, genome.fa.fai]
    versions            = ch_versions                  // channel: [versions.yml]
}