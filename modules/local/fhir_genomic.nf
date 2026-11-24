include { asFlag } from '../../subworkflows/local/utils_params'

process FHIR_GENOMIC {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12 conda-forge::openpyxl=3.1.5"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(consensus), path(serotype_json), path(genotype_json), path(coverage_json)

    output:
    tuple val(meta), path("${meta.id}.genomic.fhir.json"), emit: json
    path "versions.yml",                                   emit: versions

    script:
    def args = []
    if (asFlag(params.fhir_include_sequence)) args << "--include_sequence"
    if (params.primer_bed)                    args << "--primer_trimmed"
    if (params.performer_org_id)              args << "--performer_org_id '${params.performer_org_id}'"
    if (params.performer_org_name)            args << "--performer_org_name '${params.performer_org_name}'"
    if (params.identifier_system_base)        args << "--identifier_system_base '${params.identifier_system_base}'"
    """
    denvfhir_genomic.py \\
        --consensus ${consensus} \\
        --sample_id '${meta.id}' \\
        --serotype_json ${serotype_json} \\
        --genotype_json ${genotype_json} \\
        --coverage_json ${coverage_json} \\
        --output ${meta.id}.genomic.fhir.json \\
        --fhir_profile ${params.fhir_profile} \\
        --pipeline_version '${workflow.manifest.version}' \\
        --consensus_min_depth ${params.consensus_min_depth} \\
        --consensus_min_af ${params.consensus_min_af} \\
        --variant_min_quality ${params.variant_min_quality} \\
        --min_completeness ${params.min_completeness} \\
        ${args.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        fhir_profile: ${params.fhir_profile}
    END_VERSIONS
    """

    stub:
    """
    echo '{"resourceType": "Bundle", "type": "collection", "entry": [{"resource": {"resourceType": "DiagnosticReport", "id": "${meta.id}-genomic-report"}}]}' > ${meta.id}.genomic.fhir.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
