process BOWTIE2_BUILD_INDEX {
    tag "$meta.id"
    label 'process_high'

    // Smart caching: persist indexes across runs and pipelines in genome-specific directories
    storeDir "${params.reference_cache_dir ?: "${baseDir}/references"}/bowtie2_indexes/${fasta.baseName.replaceAll(/\.fa(sta)?$/, '')}"
    cache 'lenient'  // Use cached results even if process changes slightly

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bowtie2:2.5.4--he96a11b_6':
        'quay.io/biocontainers/bowtie2:2.5.4--he96a11b_6' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}*.bt2"), emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = meta.id
    """
    # Build bowtie2 index with organism-specific naming
    bowtie2-build \\
        --threads ${task.cpus} \\
        ${fasta} \\
        ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -n1 | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def prefix = meta.id
    """
    touch ${prefix}.1.bt2
    touch ${prefix}.2.bt2
    touch ${prefix}.3.bt2
    touch ${prefix}.4.bt2
    touch ${prefix}.rev.1.bt2
    touch ${prefix}.rev.2.bt2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(echo "2.5.4")
    END_VERSIONS
    """
}