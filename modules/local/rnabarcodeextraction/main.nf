// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process RNABARCODEEXTRACTION {
    tag "$meta.id"
    label 'process_medium'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/astral-sh/uv:python3.12-bookworm':
        'ghcr.io/astral-sh/uv:python3.12-bookworm' }"

    // Set UV cache directory to a writable location
    containerOptions = workflow.containerEngine == 'singularity' ?
        "--env UV_CACHE_DIR=\$PWD/.uv_cache" :
        "--env UV_CACHE_DIR=/tmp/.uv_cache"

    input:
    tuple val(meta), path(reads), path(config)

    output:
    tuple val(meta), path("cutadapt/*_barcoded_R{1,2}.cutadapt.fastq"), emit: barcoded_reads
    tuple val(meta), path("cutadapt/*_unbarcoded_R{1,2}.cutadapt.fastq"), emit: unbarcoded_reads
    tuple val(meta), path("cutadapt/*.cutadapt.json"), emit: json_report
    tuple val(meta), path("cutadapt/*.cutadapt.txt"), emit: txt_report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read1 = reads instanceof List ? reads[0] : reads[0]
    def read2 = reads instanceof List ? reads[1] : reads[1]

    // Check if files are gzipped and prepare file paths
    def read1_uncompressed = read1.toString().endsWith('.gz') ? read1.toString().replaceAll(/\.gz$/, '') : read1
    def read2_uncompressed = read2.toString().endsWith('.gz') ? read2.toString().replaceAll(/\.gz$/, '') : read2

    """
    # Decompress files if gzipped (paired-end mode)
    if [[ ${read1} == *.gz ]]; then
        gunzip -c ${read1} > ${read1_uncompressed}
        READ1_INPUT=${read1_uncompressed}
    else
        READ1_INPUT=${read1}
    fi

    if [[ ${read2} == *.gz ]]; then
        gunzip -c ${read2} > ${read2_uncompressed}
        READ2_INPUT=${read2_uncompressed}
    else
        READ2_INPUT=${read2}
    fi

    ${moduleDir}/templates/portable_cutadapt_wrapper.py \\
        \$READ1_INPUT \\
        \$READ2_INPUT \\
        --config-file ${config} \\
        --output-dir cutadapt \\
        --output-name-override ${prefix} \\
        --threads $task.cpus \\
        $args

    # Get version information from the wrapper
    ${moduleDir}/templates/portable_cutadapt_wrapper.py --version > version_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portable_cutadapt_wrapper: \$(grep "portable_cutadapt_wrapper.py version:" version_info.txt | cut -d' ' -f3)
        cutadapt: \$(grep "cutadapt:" version_info.txt | cut -d' ' -f2)
        python: \$(grep "Python:" version_info.txt | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p cutadapt
    touch cutadapt/${prefix}_barcoded_R1.cutadapt.fastq
    touch cutadapt/${prefix}_barcoded_R2.cutadapt.fastq
    touch cutadapt/${prefix}_unbarcoded_R1.cutadapt.fastq
    touch cutadapt/${prefix}_unbarcoded_R2.cutadapt.fastq
    touch cutadapt/${prefix}.cutadapt.json
    touch cutadapt/${prefix}.cutadapt.txt

    # Get real version information even in stub mode
    ${moduleDir}/templates/portable_cutadapt_wrapper.py --version > version_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portable_cutadapt_wrapper: \$(grep "portable_cutadapt_wrapper.py version:" version_info.txt | cut -d' ' -f3)
        cutadapt: \$(grep "cutadapt:" version_info.txt | cut -d' ' -f2)
        python: \$(grep "Python:" version_info.txt | cut -d' ' -f2)
    END_VERSIONS
    """
}
