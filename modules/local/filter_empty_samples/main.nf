process FILTER_EMPTY_SAMPLES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(reads)
    val(min_reads)

    output:
    tuple val(meta), path(reads), env(read_count), emit: samples_with_counts
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def is_paired = reads instanceof List && reads.size() == 2
    def read1 = is_paired ? reads[0] : reads
    """
    # Count reads in FASTQ file(s) (divide line count by 4)
    if [[ ${read1} == *.gz ]]; then
        read_count=\$(zcat ${read1} | wc -l | awk '{print int(\$1/4)}')
    else
        read_count=\$(wc -l < ${read1} | awk '{print int(\$1/4)}')
    fi

    echo "Sample: ${meta.id}, Reads: \${read_count}, Threshold: ${min_reads}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    read_count=1000

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n1 | sed 's/.*version //; s/ .*//' || echo "5.0")
    END_VERSIONS
    """
}
