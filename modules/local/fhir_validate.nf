include { asFlag } from '../../subworkflows/local/utils_params'

process FHIR_VALIDATE {
    label 'process_medium'
    conda "conda-forge::openjdk=21"
    container "eclipse-temurin:21-jre"

    input:
    path fhir_files, stageAs: 'input/*'
    path validator_jar
    path codesystems

    output:
    tuple path("input/*.json", includeInputs: true), path("validation.json"), emit: result
    path "validator.log",                                                    emit: log
    path "versions.yml",                                                     emit: versions

    script:
    def mem = (task.memory.toMega() * 0.8).intValue()
    """
    fhir_home="\${HOME:-}"
    if [ -z "\$fhir_home" ] || [ ! -w "\$fhir_home" ]; then
        fhir_home="\$PWD/fhir_home"
        mkdir -p "\$fhir_home"
    fi

    java -Xmx${mem}M -Dfile.encoding=UTF-8 -Duser.home="\$fhir_home" -jar ${validator_jar} \\
        -version 4.0.1 \\
        -ig ${params.fhir_ig} \\
        -ig ${codesystems} \\
        -tx ${params.fhir_tx_server} \\
        -output validation.json \\
        input/*.json \\
        > validator.log 2>&1 || true

    if [ ! -s validation.json ]; then
        echo "ERROR: FHIR validator produced no output; see validator.log" >&2
        tail -n 40 validator.log >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        java: \$(java -version 2>&1 | head -n 1 | sed 's/.*"\\(.*\\)".*/\\1/')
        fhir_ig: ${params.fhir_ig}
    END_VERSIONS
    """

    stub:
    """
    echo '{"resourceType": "Bundle", "entry": []}' > validation.json
    touch validator.log
    printf '"%s":\\n    java: stub\\n' "${task.process}" > versions.yml
    """
}

process FHIR_VALIDATION_SUMMARY {
    label 'process_single'
    conda "conda-forge::python=3.12"
    container "quay.io/biocontainers/python:3.12.12"

    input:
    tuple path(fhir_files, stageAs: 'input/*'), path(validation_json)

    output:
    path "validation_summary.tsv",            emit: summary
    path "passed/*.json",                     emit: passed, optional: true
    path "versions.yml",                      emit: versions

    script:
    def fail_default = (params.fhir_profile == 'generic')
    def fail_flag = (params.fail_on_validation_error != null) ? asFlag(params.fail_on_validation_error) : fail_default
    def fail_arg = fail_flag ? '--fail_on_error' : ''
    def offline_arg = (params.fhir_tx_server == 'n/a') ? '--offline' : ''
    """
    summarize_validation.py \\
        --validation ${validation_json} \\
        --files input/*.json \\
        ${fail_arg} \\
        ${offline_arg}

    mkdir -p passed
    while IFS= read -r fname; do
        [ -z "\$fname" ] && continue
        cp -L "input/\$fname" "passed/\$fname"
    done < passed.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    printf 'file\\terrors\\twarnings\\tinformation\\n' > validation_summary.tsv
    mkdir -p passed
    for f in input/*.json; do
        printf '%s\\t0\\t0\\t0\\n' "\$(basename "\$f")" >> validation_summary.tsv
        cp -L "\$f" "passed/\$(basename "\$f")"
    done
    printf '"%s":\\n    python: stub\\n' "${task.process}" > versions.yml
    """
}
