process PREPARE_CONSENSUS {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.consensus.fasta"), emit: consensus
    tuple val(meta), path("${meta.id}.coverage.json"),   emit: coverage
    path "versions.yml",                                 emit: versions

    script:
    """
    summarize_coverage.py \\
        --mode fasta \\
        --sample_id '${meta.id}' \\
        --input_fasta ${fasta} \\
        --output_fasta ${meta.id}.consensus.fasta \\
        --output_json ${meta.id}.coverage.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf '>${meta.id}\\nACGTACGTACGT\\n' > ${meta.id}.consensus.fasta
    echo '{"sample_id":"${meta.id}","assessed":false,"status":"submitted_consensus","completeness":1.0}' > ${meta.id}.coverage.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
