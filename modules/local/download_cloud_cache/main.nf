process DOWNLOAD_CLOUD_CACHE {
    tag "download_cloud_cache"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container 'public.ecr.aws/aws-cli/aws-cli:latest'

    // Cache reference files locally
    publishDir "${params.reference_cache_dir}/fasta", mode: params.publish_dir_mode, pattern: "*.fa"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "*.gtf"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "star_indexes/**"
    publishDir "${params.reference_cache_dir}", mode: params.publish_dir_mode, pattern: "biscuit_indexes/**"

    input:
    val genome_filename
    val gtf_filename
    val genome_id

    output:
    path "${genome_filename}", emit: genome_fasta, optional: true
    path "${gtf_filename}", emit: gtf, optional: true
    path "star_indexes/${genome_id}", emit: star_index, optional: true
    path "biscuit_indexes/${genome_id}", emit: biscuit_index, optional: true
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

    # CRITICAL: Check authentication first - fail loudly if not authenticated
    echo "🔐 Verifying cloud storage authentication..."
    if ! aws s3 ls "${params.cloud_reference_cache}/" >/dev/null 2>&1; then
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo "❌ FATAL ERROR: Cannot access cloud storage"
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        echo ""
        echo "Cloud storage fallback is enabled but authentication failed."
        echo "This would cause the pipeline to waste HOURS rebuilding indexes"
        echo "instead of downloading them in minutes."
        echo ""
        echo "⚠️  REFUSING TO CONTINUE - Please authenticate first:"
        echo ""
        echo "    aws configure"
        echo ""
        echo "Or disable cloud cache with:"
        echo ""
        echo "    --cloud_reference_cache null"
        echo ""
        echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
        exit 1
    fi
    echo "      ✓ Authentication successful"
    echo ""

    # Download genome FASTA if not already cached locally
    if [ ! -f "${params.reference_cache_dir}/fasta/${genome_filename}" ]; then
        echo "[1/4] Checking for genome FASTA: ${genome_filename}"
        if aws s3 ls "${params.cloud_reference_cache}/fasta/${genome_filename}" >/dev/null 2>&1; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/fasta/"
            echo "      ⬇  Downloading..."
            aws s3 cp "${params.cloud_reference_cache}/fasta/${genome_filename}" .
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded genome FASTA from cloud cache"
            else
                echo "      ⚠️  Download failed, will use original source"
            fi
        else
            echo "      ✗ Not found in cloud cache, will use original source"
        fi
    else
        echo "[1/4] Genome FASTA already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/fasta/${genome_filename}" .
    fi
    echo ""

    # Download GTF if not already cached locally
    if [ ! -f "${params.reference_cache_dir}/${gtf_filename}" ]; then
        echo "[2/4] Checking for GTF annotation: ${gtf_filename}"
        if aws s3 ls "${params.cloud_reference_cache}/${gtf_filename}" >/dev/null 2>&1; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/"
            echo "      ⬇  Downloading..."
            aws s3 cp "${params.cloud_reference_cache}/${gtf_filename}" .
            if [ \$? -eq 0 ]; then
                echo "      ✅ Successfully downloaded GTF annotation from cloud cache"
            else
                echo "      ⚠️  Download failed, will use original source"
            fi
        else
            echo "      ✗ Not found in cloud cache, will use original source"
        fi
    else
        echo "[2/4] GTF annotation already cached locally, skipping"
        ln -s "${params.reference_cache_dir}/${gtf_filename}" .
    fi
    echo ""

    # Download STAR index if not already cached locally
    if [ ! -d "${params.reference_cache_dir}/star_indexes/${genome_id}" ]; then
        echo "[3/4] Checking for STAR index: ${genome_id}"
        if aws s3 ls "${params.cloud_reference_cache}/star_indexes/${genome_id}/Genome" >/dev/null 2>&1; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/star_indexes/"
            echo "      ⬇  Downloading (~30GB, this may take a few minutes)..."
            mkdir -p star_indexes
            aws s3 cp --recursive "${params.cloud_reference_cache}/star_indexes/${genome_id}" "star_indexes/${genome_id}/"
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
        echo "[3/4] STAR index already cached locally, skipping"
    fi
    echo ""

    # Download Biscuit index if not already cached locally
    if [ ! -d "${params.reference_cache_dir}/biscuit_indexes/${genome_id}" ]; then
        echo "[4/4] Checking for Biscuit index: ${genome_id}"
        if aws s3 ls "${params.cloud_reference_cache}/biscuit_indexes/${genome_id}/${genome_filename}.dau.bwt" >/dev/null 2>&1; then
            echo "      ✓ Found in cloud cache at: ${params.cloud_reference_cache}/biscuit_indexes/"
            echo "      ⬇  Downloading (~13GB, this may take a few minutes)..."
            mkdir -p biscuit_indexes
            aws s3 cp --recursive "${params.cloud_reference_cache}/biscuit_indexes/${genome_id}" "biscuit_indexes/${genome_id}/"
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
        echo "[4/4] Biscuit index already cached locally, skipping"
    fi
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(aws --version 2>&1 | cut -d' ' -f1 | cut -d'/' -f2)
    END_VERSIONS
    """

    stub:
    """
    # Create mock files for testing
    touch ${genome_filename}
    touch ${gtf_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(echo "2.15.0")
    END_VERSIONS
    """
}
