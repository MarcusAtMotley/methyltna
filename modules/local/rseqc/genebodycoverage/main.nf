process RSEQC_GENEBODYCOVERAGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rseqc:5.0.3--py39hf95cd2a_0' :
        'biocontainers/rseqc:5.0.3--py39hf95cd2a_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
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
    # Run geneBody_coverage.py with error handling for low-coverage samples
    # If the tool fails (e.g., insufficient coverage), create minimal output to allow pipeline continuation
    geneBody_coverage.py \\
        -i $bam \\
        -r $bed \\
        -o ${prefix} \\
        $args || {
        echo "WARNING: geneBody_coverage.py failed for ${prefix} (likely insufficient mapped reads)" >&2
        echo "Creating minimal output to allow pipeline continuation" >&2
        # Remove any partially created files before creating fallback output (handles bash noclobber mode)
        rm -f ${prefix}.geneBodyCoverage.txt ${prefix}.geneBodyCoverage.curves.pdf
        echo "# Gene body coverage analysis failed - insufficient mapped reads or coverage" > ${prefix}.geneBodyCoverage.txt
        echo "# This sample will be excluded from gene body coverage visualization in MultiQC" >> ${prefix}.geneBodyCoverage.txt

        # Create versions file for failed samples too
        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rseqc: \$(geneBody_coverage.py --version | sed -e "s/geneBody_coverage.py //g")
    END_VERSIONS

        exit 0
    }

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
