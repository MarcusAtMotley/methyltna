process DOWNLOAD_CLOUD_CACHE {
    tag "download_cloud_cache"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'biocontainers/gcloud:498.0.0':
        'quay.io/biocontainers/gcloud:498.0.0' }"

    // Cache reference files locally
    publishDir "${params.reference_cache_dir}/fasta", mode: params.publish_dir_mode, pattern: "*.fa"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "*.gtf"

    input:
    val genome_filename
    val gtf_filename

    output:
    path "${genome_filename}", emit: genome_fasta, optional: true
    path "${gtf_filename}", emit: gtf, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Try to download from cloud reference cache
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "📦 Checking cloud reference cache..."
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""

    # Download genome FASTA if not already cached locally
    if [ ! -f "${params.reference_cache_dir}/fasta/${genome_filename}" ]; then
        echo "[1/2] Checking for genome FASTA: ${genome_filename}"
        if gsutil -q stat "${params.cloud_reference_cache}/fasta/${genome_filename}" 2>/dev/null; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/fasta/"
            echo "      ⬇  Downloading..."
            gsutil -m cp "${params.cloud_reference_cache}/fasta/${genome_filename}" .
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded genome FASTA from cloud cache"
            else
                echo "      ⚠️  Download failed, will use original source"
            fi
        else
            echo "      ✗ Not found in cloud cache, will use original source"
        fi
    else
        echo "[1/2] Genome FASTA already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/fasta/${genome_filename}" .
    fi
    echo ""

    # Download GTF if not already cached locally
    if [ ! -f "${params.reference_cache_dir}/${gtf_filename}" ]; then
        echo "[2/2] Checking for GTF annotation: ${gtf_filename}"
        if gsutil -q stat "${params.cloud_reference_cache}/${gtf_filename}" 2>/dev/null; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/"
            echo "      ⬇  Downloading..."
            gsutil -m cp "${params.cloud_reference_cache}/${gtf_filename}" .
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded GTF annotation from cloud cache"
            else
                echo "      ⚠️  Download failed, will use original source"
            fi
        else
            echo "      ✗ Not found in cloud cache, will use original source"
        fi
    else
        echo "[2/2] GTF annotation already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/${gtf_filename}" .
    fi
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gsutil: \$(gsutil version | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
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
