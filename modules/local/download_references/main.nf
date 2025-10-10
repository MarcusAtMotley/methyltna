process DOWNLOAD_REFERENCES {
    tag "download_references"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/gcloud:498.0.0':
        'quay.io/biocontainers/gcloud:498.0.0' }"

    // Cache reference files to avoid re-downloading
    publishDir "${params.reference_cache_dir}/fasta", mode: params.publish_dir_mode, pattern: "*.fa"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "*.gtf"

    output:
    path "${params.genome_fasta.tokenize('/').last()}", emit: genome_fasta
    path "${params.annotation_gtf.tokenize('/').last()}", emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def genome_filename = params.genome_fasta.tokenize('/').last()
    def gtf_filename = params.annotation_gtf.tokenize('/').last()
    """
    # Function to provide helpful error messages for authentication failures
    download_with_error_handling() {
        local url="\$1"
        local output="\$2"
        local description="\$3"

        echo "Downloading \$description from \$url..."

        if ! gsutil cp "\$url" "\$output" 2>/dev/null; then
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "❌ DOWNLOAD ERROR: Failed to download \$description"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo ""
            echo "🔧 MANUAL DOWNLOAD REQUIRED:"
            echo "The pipeline cannot download the required reference files. This could be due to:"
            echo "  • Authentication issues with Google Cloud Storage"
            echo "  • Network connectivity problems"
            echo "  • Container access issues"
            echo ""
            echo "📥 Please download this file manually:"
            echo "   Source: \$url"
            echo "   Save to: ${params.reference_cache_dir}/fasta/\$output"
            echo ""
            echo "💡 Commands to run:"
            echo "   mkdir -p ${params.reference_cache_dir}/fasta"
            echo "   gsutil cp \$url ${params.reference_cache_dir}/fasta/\$output"
            echo ""
            echo "   OR if you have the file locally:"
            echo "   cp /path/to/your/\$output ${params.reference_cache_dir}/fasta/\$output"
            echo ""
            echo "🔄 Then rerun the pipeline with: nextflow run . -resume <other_params>"
            echo ""
            echo "⚙️  Pipeline Parameters:"
            echo "   --force_redownload_references   Force re-download even if files are cached"
            echo "   --force_rebuild_indexes         Force rebuild even if indexes are cached"
            echo ""
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            return 1
        fi

        echo "✅ Successfully downloaded \$description"
        return 0
    }

    # Download genome reference with error handling
    download_with_error_handling "${params.genome_fasta}" "${genome_filename}" "genome reference (${genome_filename})"
    if [ \$? -ne 0 ]; then exit 1; fi

    # Download GTF file with error handling
    download_with_error_handling "${params.annotation_gtf}" "${gtf_filename}" "GTF annotation (${gtf_filename})"
    if [ \$? -ne 0 ]; then exit 1; fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsutil: \$(gsutil version | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def genome_filename = params.genome_fasta.tokenize('/').last()
    def gtf_filename = params.annotation_gtf.tokenize('/').last()
    """
    # Create mock files for testing
    touch ${genome_filename}
    touch ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsutil: \$(echo "498.0.0")
    END_VERSIONS
    """
}