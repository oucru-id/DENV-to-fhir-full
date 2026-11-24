process NANOPLOT {
    tag "${meta.id}"
    label 'process_medium'
    conda "bioconda::nanoplot=1.48.0"
    container "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_NanoStats.txt"),        emit: stats
    tuple val(meta), path("${meta.id}_NanoPlot-report.html"), emit: html
    path "versions.yml",                                      emit: versions

    script:
    """
    NanoPlot \\
        --fastq ${reads} \\
        --threads ${task.cpus} \\
        --outdir . \\
        --prefix ${meta.id}_ \\
        --no_static \\
        --N50

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(NanoPlot --version 2>&1 | sed 's/NanoPlot //')
    END_VERSIONS
    """

    stub:
    """
    printf 'Total bases:\\t100,000\\n' > ${meta.id}_NanoStats.txt
    touch ${meta.id}_NanoPlot-report.html
    printf '"%s":\\n    nanoplot: stub\\n' "${task.process}" > versions.yml
    """
}
