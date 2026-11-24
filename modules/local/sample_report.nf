process SAMPLE_REPORT {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(consensus), path(serotype_json), path(genotype_json), path(coverage_json)

    output:
    tuple val(meta), path("${meta.id}.summary_report.txt"), emit: report
    path "versions.yml",                                    emit: versions

    script:
    """
    generate_sample_report.py \\
        --sample_id '${meta.id}' \\
        --consensus ${consensus} \\
        --serotype_json ${serotype_json} \\
        --genotype_json ${genotype_json} \\
        --coverage_json ${coverage_json} \\
        --output ${meta.id}.summary_report.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'DENVtoFHIR Summary Report: ${meta.id}\\n' > ${meta.id}.summary_report.txt
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
