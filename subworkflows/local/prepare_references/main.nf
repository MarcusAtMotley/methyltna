#!/usr/bin/env nextflow

//
// PREPARE_REFERENCES: Smart reference download and indexing with caching
//

include { DOWNLOAD_REFERENCES } from '../../../modules/local/download_references/main'
include { DOWNLOAD_CLOUD_CACHE } from '../../../modules/local/download_cloud_cache/main'
include { STAR_GENOMEGENERATE } from '../../../modules/nf-core/star/genomegenerate/main'
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

    // Check for locally cached files and indexes
    def genome_cached = file("${params.reference_cache_dir}/fasta/${genome_filename}").exists()
    def gtf_cached = file("${params.reference_cache_dir}/${gtf_filename}").exists()
    def star_index_cached = file("${params.reference_cache_dir}/star_indexes/${genome_id}").exists() &&
                            file("${params.reference_cache_dir}/star_indexes/${genome_id}/Genome").exists()
    def biscuit_index_cached = file("${params.reference_cache_dir}/biscuit_indexes/${genome_id}").exists() &&
                               file("${params.reference_cache_dir}/biscuit_indexes/${genome_id}/${genome_filename}.dau.bwt").exists()

    // Determine if we need to download based on cache status, provided paths, and force flags
    def need_download = false
    def need_cloud_download = false

    // Check if we need FASTA/GTF from cloud or original source
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

    // Also check if any indexes are missing and cloud cache is available
    if (!params.force_rebuild_indexes && params.cloud_reference_cache) {
        if (!star_index_cached || !biscuit_index_cached) {
            need_cloud_download = true
        }
    }

    // Placeholder channels for indexes that will be conditionally filled
    ch_star_index_from_cloud = Channel.empty()
    ch_biscuit_index_from_cloud = Channel.empty()

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

        // Get index channels from cloud download (may be empty if not found)
        ch_star_index_from_cloud = DOWNLOAD_CLOUD_CACHE.out.star_index
            .map { dir -> [[id:'genome'], dir] }
        ch_biscuit_index_from_cloud = DOWNLOAD_CLOUD_CACHE.out.biscuit_index
            .map { dir -> [[id:'genome'], dir] }

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

    // Build or use STAR genome index - priority: user-provided > cloud > local cache > build
    // Use pure channel operations to ensure proper dependency ordering
    def star_index_dir = "${params.reference_cache_dir}/star_indexes/${genome_id}"

    if (params.star_index != null) {
        // Priority 1: User provided explicit STAR index path
        ch_star_index = Channel.fromPath(params.star_index)
            .map { dir -> [[id:'genome'], dir] }
    } else {
        // Priority 2-4: Cloud download > Local cache > Build from scratch
        // Use .ifEmpty() to create runtime dependency chain
        ch_star_index = ch_star_index_from_cloud
            .ifEmpty {
                // Cloud didn't provide - check local cache at runtime
                def star_cached = file("${star_index_dir}").exists() &&
                                  file("${star_index_dir}/Genome").exists()

                if (star_cached && !params.force_rebuild_indexes) {
                    // Priority 3: Use local cache
                    Channel.fromPath(star_index_dir).map { dir -> [[id:'genome'], dir] }
                } else {
                    // Priority 4: Build from scratch
                    def ch_genome_with_id = ch_genome_fasta.map { meta, fasta ->
                        [meta + [genome_id: genome_id], fasta]
                    }
                    def ch_gtf_with_id = ch_annotation_gtf.map { meta, gtf ->
                        [meta + [genome_id: genome_id], gtf]
                    }
                    STAR_GENOMEGENERATE(ch_genome_with_id, ch_gtf_with_id)
                    ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
                    STAR_GENOMEGENERATE.out.index
                }
            }
    }


    // Build or use Biscuit genome index - priority: user-provided > cloud > local cache > build
    // Use pure channel operations to ensure proper dependency ordering
    def biscuit_index_dir = "${params.reference_cache_dir}/biscuit_indexes/${genome_id}"

    if (params.biscuit_index != null) {
        // Priority 1: User provided explicit Biscuit index path
        ch_biscuit_index = Channel.fromPath(params.biscuit_index)
            .map { dir -> [[id:'genome'], dir] }
    } else {
        // Priority 2-4: Cloud download > Local cache > Build from scratch
        // Use .ifEmpty() to create runtime dependency chain
        ch_biscuit_index = ch_biscuit_index_from_cloud
            .ifEmpty {
                // Cloud didn't provide - check local cache at runtime
                def biscuit_cached = file("${biscuit_index_dir}").exists() &&
                                     file("${biscuit_index_dir}/${genome_filename}.dau.bwt").exists()

                if (biscuit_cached && !params.force_rebuild_indexes) {
                    // Priority 3: Use local cache
                    Channel.fromPath(biscuit_index_dir).map { dir -> [[id:'genome'], dir] }
                } else {
                    // Priority 4: Build from scratch
                    def ch_genome_with_id = ch_genome_fasta.map { meta, fasta ->
                        [meta + [genome_id: genome_id], fasta]
                    }
                    BISCUIT_INDEX(ch_genome_with_id)
                    ch_versions = ch_versions.mix(BISCUIT_INDEX.out.versions)
                    BISCUIT_INDEX.out.index
                }
            }
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
    biscuit_index       = ch_biscuit_index             // channel: [meta, biscuit_index_dir]
    genome_fai          = ch_genome_fai                // channel: [meta, genome.fa.fai]
    versions            = ch_versions                  // channel: [versions.yml]
}