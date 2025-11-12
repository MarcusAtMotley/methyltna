process SUBSAMPLE_BAM_FOR_QC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    val(max_reads)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Count mapped reads (exclude unmapped and secondary alignments)
    # Flag 260 = unmapped (4) + not primary alignment (256)
    total_reads=\$(samtools view -c -F 260 ${bam})

    echo "Sample: ${prefix}"
    echo "  Total mapped reads: \${total_reads}"

    if [ \${total_reads} -gt ${max_reads} ]; then
        # Calculate subsampling fraction
        fraction=\$(awk "BEGIN {printf \\"%.6f\\", ${max_reads}/\${total_reads}}")
        echo "  Subsampling to ${max_reads} reads (fraction: \${fraction})"

        # Subsample BAM file with random seed for reproducibility
        # Strip leading "0." from fraction for samtools -s format (seed.fraction)
        samtools view -s 42.\${fraction#0.} -b ${args} ${bam} > ${prefix}_subsampled.bam

        # Index subsampled BAM
        samtools index ${prefix}_subsampled.bam
    else
        # Use full BAM (below threshold)
        echo "  Using all reads (below ${max_reads} threshold)"
        ln -s ${bam} ${prefix}_subsampled.bam
        ln -s ${bai} ${prefix}_subsampled.bam.bai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_subsampled.bam
    touch ${prefix}_subsampled.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' || echo "1.21")
    END_VERSIONS
    """
}
