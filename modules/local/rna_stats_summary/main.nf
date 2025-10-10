process RNA_STATS_SUMMARY {
    tag "RNA barcode statistics"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/astral-sh/uv:python3.12-bookworm':
        'ghcr.io/astral-sh/uv:python3.12-bookworm' }"

    // Set UV cache directory to a writable location
    containerOptions = workflow.containerEngine == 'singularity' ?
        "--env UV_CACHE_DIR=\$PWD/.uv_cache" :
        "--env UV_CACHE_DIR=/tmp/.uv_cache"

    input:
    path(json_files)

    output:
    path("rna_barcode_stats.tsv"), emit: tsv
    path("rna_barcode_stats_mqc.txt"), emit: mqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    ${projectDir}/bin/extract_rna_stats.py ${json_files.join(' ')} -o rna_barcode_stats.tsv $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        uv: \$(uv --version 2>&1 | sed 's/uv //' || echo "Unknown")
        python: \$(python3 --version 2>&1 | sed 's/Python //' || echo "Unknown")
        extract_rna_stats: 1.0.0
    END_VERSIONS
    """

    stub:
    """
    touch rna_barcode_stats.tsv
    touch rna_barcode_stats_mqc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //' || echo "3.12.0")
        extract_rna_stats: 1.0.0
    END_VERSIONS
    """
}