process RSEQC_GENEBODYCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam)
    path  bed

    output:
    tuple val(meta), path("*.geneBodyCoverage.txt")       , emit: txt
    tuple val(meta), path("*.geneBodyCoverage.curves.pdf"), emit: pdf, optional: true
    path  "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    geneBody_coverage.py \\
        -i $bam \\
        -r $bed \\
        -o ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(geneBody_coverage.py --version | sed -e "s/geneBody_coverage.py //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.geneBodyCoverage.txt
    touch ${prefix}.geneBodyCoverage.curves.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: 5.0.3
    END_VERSIONS
    """
}
