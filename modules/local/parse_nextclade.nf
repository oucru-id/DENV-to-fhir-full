process PARSE_NEXTCLADE {
    tag "${meta.id}"
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple val(meta), path(consensus), path(serotype_json), path(coverage_json), path(nextclade_csv)

    output:
    tuple val(meta), path("${meta.id}.genotype.json"), emit: json
    path "versions.yml",                               emit: versions

    script:
    def csv_arg = nextclade_csv.name != 'NO_FILE' ? "--nextclade_csv ${nextclade_csv}" : ""
    """
    parse_nextclade.py \\
        --sample_id '${meta.id}' \\
        --consensus ${consensus} \\
        --serotype_json ${serotype_json} \\
        --coverage_json ${coverage_json} \\
        --min_completeness ${params.min_completeness} \\
        --output ${meta.id}.genotype.json \\
        ${csv_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample_id":"${meta.id}","serotype":"DENV-1","genotype":"Genotype V","major_lineage":"Lineage D","minor_lineage":"Sub-lineage 1.1","clade":"1V_D.1.1","qc_status":"PASS","qc_reasons":[],"confidence":"high","aa_substitutions":[],"aa_deletions":[],"aa_insertions":[]}' > ${meta.id}.genotype.json
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
