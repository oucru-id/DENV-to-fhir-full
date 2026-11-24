process NEXTCLADE_DATASET_GET {
    tag "${dataset_name}"
    label 'process_low'
    conda "bioconda::nextclade=3.29.1"
    container "quay.io/biocontainers/nextclade:3.29.1--h9ee0642_0"

    input:
    val dataset_name

    output:
    tuple val(dataset_name), path("nextclade_dataset"), emit: dataset
    path "versions.yml",                                emit: versions

    script:
    def tag_arg = params.nextclade_dataset_tag ? "--tag '${params.nextclade_dataset_tag}'" : ""
    """
    nextclade dataset get \\
        --name '${dataset_name}' \\
        ${tag_arg} \\
        --output-dir nextclade_dataset

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nextclade: \$(nextclade --version | sed 's/nextclade //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p nextclade_dataset
    echo '{"name":"${dataset_name}"}' > nextclade_dataset/pathogen.json
    printf '"%s":\\n    nextclade: stub\\n' "${task.process}" > versions.yml
    """
}
