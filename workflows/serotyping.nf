#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_MINIMAP2_INDEX {
    input:
    path(references)
    
    output:
    path("denv_ref.mmi"),                       emit: index
    path("all_denv_references_tagged.fasta"),    emit: tagged_fasta
    
    script:
    """
    if [ \$(ls -1 ${references} 2>/dev/null | wc -l) -eq 0 ]; then
        echo "No reference files found"
        exit 1
    fi

    cat ${references} > all_denv_references_raw.fasta

    if [ ! -s all_denv_references_raw.fasta ]; then
        echo "Combined reference file is empty"
        exit 1
    fi

    if ! grep -q "^>" all_denv_references_raw.fasta; then
        echo "No FASTA headers found in reference file"
        exit 1
    fi

    SEQ_COUNT=\$(grep -c "^>" all_denv_references_raw.fasta)
    if [ \$SEQ_COUNT -lt 4 ]; then
        echo "WARNING: Only \$SEQ_COUNT reference sequences found"
    fi

    sed \\
        's/>NC_001477\\.1/>DENV1|NC_001477.1/g;
         s/>NC_001474\\.2/>DENV2|NC_001474.2/g;
         s/>NC_001475\\.2/>DENV3|NC_001475.2/g;
         s/>NC_002640\\.1/>DENV4|NC_002640.1/g;
         s/>JF262780\\.1/>DENV4|JF262780.1/g;
         s/>EF105380\\.1/>DENV2_SYLVATIC|EF105380.1/g' \\
        all_denv_references_raw.fasta > all_denv_references_tagged.fasta

    minimap2 -d denv_ref.mmi all_denv_references_tagged.fasta

    if [ ! -s denv_ref.mmi ]; then
        echo "minimap2 index creation failed"
        exit 1
    fi
    """
}

process SEROTYPE_FROM_READS {
    publishDir "${params.results_dir}/serotyping", mode: 'copy', pattern: "*.serotype.json"
    
    input:
    tuple val(sample_id), path(reads)
    path(minimap2_index)
    
    output:
    tuple val(sample_id), path("${sample_id}.serotype.json"),   emit: serotype_info
    tuple val(sample_id), path("${sample_id}.minimap2.paf"),    emit: minimap2_paf
    
    script:
    def is_illumina = reads instanceof List
    def read_file   = is_illumina ? reads[0] : reads
    def mm2_preset  = is_illumina ? "sr" : "map-ont"
    """
    create_fallback_json() {
        echo "{\\"sample_id\\": \\"${sample_id}\\", \\"serotype\\": \\"unknown\\", \\"confidence\\": \\"none\\", \\"note\\": \\"\$1\\"}" > ${sample_id}.serotype.json
    }

    if [ ! -s "${read_file}" ]; then
        touch ${sample_id}.minimap2.paf
        create_fallback_json "empty_input_file"
        exit 0
    fi

    if [[ "${read_file}" == *.gz ]]; then
        FIRST_LINE=\$(zcat "${read_file}" | head -n 1)
    else
        FIRST_LINE=\$(head -n 1 "${read_file}")
    fi

    if [[ ! "\$FIRST_LINE" =~ ^@ ]]; then
        touch ${sample_id}.minimap2.paf
        create_fallback_json "invalid_fastq_header"
        exit 0
    fi

    INITIAL_READS=60000
    MAX_READS=180000

    for ATTEMPT in 1 2 3; do
        CURRENT_READS=\$((INITIAL_READS * ATTEMPT * ATTEMPT))
        if [ \$CURRENT_READS -gt \$MAX_READS ]; then
            CURRENT_READS=\$MAX_READS
        fi

        LINES=\$((CURRENT_READS * 4))

        if [[ "${read_file}" == *.gz ]]; then
            zcat "${read_file}" | head -n \$LINES > subset.fastq 2>/dev/null || zcat "${read_file}" > subset.fastq
        else
            head -n \$LINES "${read_file}" > subset.fastq 2>/dev/null || cat "${read_file}" > subset.fastq
        fi

        ACTUAL_READS=\$(( \$(wc -l < subset.fastq) / 4 ))
        echo "Attempt \$ATTEMPT: extracted \$ACTUAL_READS reads"

        if [ \$ACTUAL_READS -eq 0 ]; then
            touch ${sample_id}.minimap2.paf
            create_fallback_json "no_reads_extracted"
            exit 0
        fi

        minimap2 -x ${mm2_preset} --secondary=no -t ${task.cpus} \\
            denv_ref.mmi subset.fastq > ${sample_id}.minimap2.paf || touch ${sample_id}.minimap2.paf

        python3 $baseDir/scripts/serotype_classification.py \\
            --blast_results ${sample_id}.minimap2.paf \\
            --format paf \\
            --sample_id ${sample_id} \\
            --total_reads \$ACTUAL_READS \\
            --output ${sample_id}.serotype.json

        if [ ! -f "${sample_id}.serotype.json" ]; then
            echo "ERROR: Classification failed"
            create_fallback_json "classification_script_failed"
            break
        fi

        CONFIDENCE=\$(python3 -c "
import json
try:
    with open('${sample_id}.serotype.json') as f:
        data = json.load(f)
    print(data.get('confidence', 'low'))
except Exception:
    print('low')
")

        if [[ "\$CONFIDENCE" == "high" ]] || [[ "\$CONFIDENCE" == "medium" ]] || [[ \$CURRENT_READS -ge \$MAX_READS ]]; then
            break
        fi

        echo "Low confidence (\$CONFIDENCE), retrying with more reads"
    done

    if [ ! -f "${sample_id}.serotype.json" ]; then
        create_fallback_json "process_failed"
    fi

    cat ${sample_id}.serotype.json
    """
}

workflow SEROTYPING {
    take:
    reads_ch
    
    main:
    references     = Channel.fromPath("${params.reference_dir}/*.fasta").collect()
    minimap2_index = CREATE_MINIMAP2_INDEX(references)
    
    SEROTYPE_FROM_READS(reads_ch, minimap2_index.index)
    
    emit:
    serotype_info  = SEROTYPE_FROM_READS.out.serotype_info
    minimap2_index = minimap2_index.index
}
