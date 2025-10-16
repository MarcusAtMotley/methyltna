process DOWNLOAD_CLOUD_CACHE {
    tag "download_cloud_cache"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'docker://google/cloud-sdk:alpine'

    // Cache reference files locally
    publishDir "${params.reference_cache_dir}/fasta", mode: params.publish_dir_mode, pattern: "*.fa"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "*.gtf"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "star_indexes"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "biscuit_indexes"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "bowtie2_indexes"

    input:
    val genome_filename
    val gtf_filename
    val genome_id

    output:
    path "${genome_filename}", emit: genome_fasta, optional: true
    path "${gtf_filename}", emit: gtf, optional: true
    path "star_indexes", emit: star_index, optional: true
    path "biscuit_indexes", emit: biscuit_index, optional: true
    path "bowtie2_indexes", emit: bowtie2_index, optional: true
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
        echo "[1/5] Checking for genome FASTA: ${genome_filename}"
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
        echo "[1/5] Genome FASTA already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/fasta/${genome_filename}" .
    fi
    echo ""

    # Download GTF if not already cached locally
    if [ ! -f "${params.reference_cache_dir}/${gtf_filename}" ]; then
        echo "[2/5] Checking for GTF annotation: ${gtf_filename}"
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
        echo "[2/5] GTF annotation already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/${gtf_filename}" .
    fi
    echo ""

    # Download STAR index if not already cached locally
    if [ ! -d "${params.reference_cache_dir}/star_indexes/${genome_id}" ]; then
        echo "[3/5] Checking for STAR index: ${genome_id}"
        if gsutil -q stat "${params.cloud_reference_cache}/star_indexes/${genome_id}/Genome" 2>/dev/null; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/star_indexes/"
            echo "      ⬇  Downloading (~30GB, this may take a few minutes)..."
            mkdir -p star_indexes
            gsutil -m cp -r "${params.cloud_reference_cache}/star_indexes/${genome_id}" star_indexes/
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded STAR index from cloud cache"
            else
                echo "      ⚠️  Download failed, will build from scratch"
                rm -rf star_indexes
            fi
        else
            echo "      ✗ Not found in cloud cache, will build from scratch"
        fi
    else
        echo "[3/5] STAR index already cached locally, skipping"
    fi
    echo ""

    # Download Biscuit index if not already cached locally
    if [ ! -d "${params.reference_cache_dir}/biscuit_indexes/${genome_id}" ]; then
        echo "[4/5] Checking for Biscuit index: ${genome_id}"
        if gsutil -q stat "${params.cloud_reference_cache}/biscuit_indexes/${genome_id}/${genome_filename}.dau.bwt" 2>/dev/null; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/biscuit_indexes/"
            echo "      ⬇  Downloading (~13GB, this may take a few minutes)..."
            mkdir -p biscuit_indexes
            gsutil -m cp -r "${params.cloud_reference_cache}/biscuit_indexes/${genome_id}" biscuit_indexes/
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded Biscuit index from cloud cache"
            else
                echo "      ⚠️  Download failed, will build from scratch"
                rm -rf biscuit_indexes
            fi
        else
            echo "      ✗ Not found in cloud cache, will build from scratch"
        fi
    else
        echo "[4/5] Biscuit index already cached locally, skipping"
    fi
    echo ""

    # Download Bowtie2 index if not already cached locally
    if [ ! -d "${params.reference_cache_dir}/bowtie2_indexes/${genome_id}" ]; then
        echo "[5/5] Checking for Bowtie2 index: ${genome_id}"
        if gsutil -q stat "${params.cloud_reference_cache}/bowtie2_indexes/${genome_id}/genome.1.bt2" 2>/dev/null; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/bowtie2_indexes/"
            echo "      ⬇  Downloading (~4GB, this may take a minute)..."
            mkdir -p bowtie2_indexes
            gsutil -m cp -r "${params.cloud_reference_cache}/bowtie2_indexes/${genome_id}" bowtie2_indexes/
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded Bowtie2 index from cloud cache"
            else
                echo "      ⚠️  Download failed, will build from scratch"
                rm -rf bowtie2_indexes
            fi
        else
            echo "      ✗ Not found in cloud cache, will build from scratch"
        fi
    else
        echo "[5/5] Bowtie2 index already cached locally, skipping"
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
