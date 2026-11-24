process MOSDEPTH {
    tag "${meta.id}"
    label 'process_low'
    conda "bioconda::mosdepth=0.3.6"
    container "quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0"

    input:
    tuple val(meta), path(bam), path(bai), path(reference)

    output:
    tuple val(meta), path("${meta.id}.per-base.bed.gz"), path("${meta.id}.mosdepth.summary.txt"), path(reference), emit: coverage
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"),                                                      emit: mqc
    path "versions.yml",                                                                                           emit: versions

    script:
    """
    mosdepth \\
        --threads ${task.cpus} \\
        ${meta.id} \\
        ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version 2>&1 | sed 's/mosdepth //')
    END_VERSIONS
    """

    stub:
    """
    printf 'NC_001477.1\\t0\\t10735\\t50\\n' | gzip > ${meta.id}.per-base.bed.gz
    printf 'chrom\\tlength\\tbases\\tmean\\tmin\\tmax\\nNC_001477.1\\t10735\\t536750\\t50.0\\t50\\t50\\n' > ${meta.id}.mosdepth.summary.txt
    printf '"%s":\\n    mosdepth: stub\\n' "${task.process}" > versions.yml
    """
}
