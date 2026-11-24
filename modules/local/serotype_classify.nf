process SEROTYPE_CLASSIFY {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(paf), val(total_reads)
    path reference_panel

    output:
    tuple val(meta), path("${meta.id}.serotype.json"),                                             emit: json
    tuple val(meta), env('DENV_SEROTYPE'), env('DENV_REF_FILE'), env('DENV_NEXTCLADE_DATASET'),    emit: assignment
    path "versions.yml",                                                                           emit: versions

    script:
    def min_reads = meta.platform == 'fasta' ? 1 : params.min_serotype_reads
    """
    serotype_classifier.py \\
        --paf ${paf} \\
        --panel ${reference_panel} \\
        --sample_id '${meta.id}' \\
        --platform '${meta.platform}' \\
        --total_reads ${total_reads} \\
        --min_reads ${min_reads} \\
        --coinfection_min_fraction ${params.coinfection_min_fraction} \\
        --output ${meta.id}.serotype.json \\
        --env_output assignment.env

    source assignment.env

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","serotype":"DENV-1","confidence":"high","selected_reference":"NC_0014771.fasta","nextclade_dataset":"community/v-gen-lab/dengue/denv1"}' > ${meta.id}.serotype.json
    DENV_SEROTYPE="DENV-1"
    DENV_REF_FILE="NC_0014771.fasta"
    DENV_NEXTCLADE_DATASET="community/v-gen-lab/dengue/denv1"
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
