#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_BLAST_DB {
    input:
    path(references)
    
    output:
    path("denv_db*"), emit: blast_db
    
    script:
    """
    if [ \$(ls -1 ${references} 2>/dev/null | wc -l) -eq 0 ]; then
        exit 1
    fi
    
    cat ${references} > all_denv_references.fasta
    
    if [ ! -s all_denv_references.fasta ]; then
        echo "Combined reference file is empty"
        exit 1
    fi
    
    if ! grep -q "^>" all_denv_references.fasta; then
        echo "No FASTA headers found in reference file"
        exit 1
    fi
    
    SEQ_COUNT=\$(grep -c "^>" all_denv_references.fasta)
    
    if [ \$SEQ_COUNT -lt 4 ]; then
        echo "WARNING: Only \$SEQ_COUNT reference sequences found"
    fi
    
    makeblastdb -in all_denv_references.fasta -dbtype nucl -out denv_db
    
    if [ ! -f denv_db.ndb ] && [ ! -f denv_db.ndb ]; then
        echo "BLAST database creation failed"
        exit 1
    fi
    """
}

process SEROTYPE_FROM_READS {
    publishDir "${params.results_dir}/serotyping", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reads)
    path(blast_db)
    
    output:
    tuple val(sample_id), path("${sample_id}.serotype.json"), emit: serotype_info
    tuple val(sample_id), path("${sample_id}.quick_blast.txt"), emit: blast_results
    
    script:
    def read_file = reads instanceof List ? reads[0] : reads
    """
    create_fallback_json() {
        echo "{\\"sample_id\\": \\"${sample_id}\\", \\"serotype\\": \\"unknown\\", \\"confidence\\": \\"none\\", \\"note\\": \\"\$1\\"}" > ${sample_id}.serotype.json
    }

    if [ ! -s "${read_file}" ]; then
        touch ${sample_id}.quick_blast.txt
        create_fallback_json "empty_input_file"
        exit 0
    fi
    
    if [[ "${read_file}" == *.gz ]]; then
        FIRST_LINE=\$(zcat ${read_file} | head -n 1)
    else
        FIRST_LINE=\$(head -n 1 ${read_file})
    fi
    
    if [[ ! "\$FIRST_LINE" =~ ^@ ]]; then
        touch ${sample_id}.quick_blast.txt
        create_fallback_json "invalid_fastq_header"
        exit 0
    fi
    
    INITIAL_READS=300000
    MAX_READS=500000  
    
    for ATTEMPT in 1 2 3; do
        CURRENT_READS=\$((INITIAL_READS * ATTEMPT * ATTEMPT))
        if [ \$CURRENT_READS -gt \$MAX_READS ]; then
            CURRENT_READS=\$MAX_READS
        fi
        
        LINES=\$((CURRENT_READS * 4))
        
        if [[ "${read_file}" == *.gz ]]; then
            zcat ${read_file} | head -n \$LINES > subset.fastq 2>/dev/null || {
                zcat ${read_file} > subset.fastq
            }
        else
            head -n \$LINES ${read_file} > subset.fastq 2>/dev/null || {
                cat ${read_file} > subset.fastq
            }
        fi
        
        ACTUAL_READS=\$(( \$(wc -l < subset.fastq) / 4 ))
        echo "Extracted \$ACTUAL_READS reads for analysis"
        
        if [ \$ACTUAL_READS -eq 0 ]; then
            touch ${sample_id}.quick_blast.txt
            create_fallback_json "no_reads_extracted"
            exit 0
        fi
        
        awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)} NR%4==2{print}' subset.fastq > subset.fasta
        
        FASTA_SEQS=\$(grep -c "^>" subset.fasta || echo 0)
        
        if [ \$FASTA_SEQS -eq 0 ]; then
            continue
        fi
        
        blastn -query subset.fasta \\
            -db denv_db \\
            -out ${sample_id}.quick_blast.txt \\
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \\
            -max_target_seqs 1 \\
            -evalue 1e-3 \\
            -word_size 11 \\
            -num_threads ${task.cpus} || touch ${sample_id}.quick_blast.txt
        
        python3 $baseDir/scripts/serotype_classification.py \\
            --blast_results ${sample_id}.quick_blast.txt \\
            --sample_id ${sample_id} \\
            --total_reads \$ACTUAL_READS \\
            --output ${sample_id}.serotype.json
        
        if [ ! -f "${sample_id}.serotype.json" ]; then
            echo "ERROR: Classification failed"
            create_fallback_json "classification_script_failed"
            break
        fi
        
        CONFIDENCE=\$(python3 -c "
import json, sys
try:
    with open('${sample_id}.serotype.json', 'r') as f:
        data = json.load(f)
    print(data.get('confidence', 'low'))
except Exception as e:
    print('low')
")
                
        if [[ "\$CONFIDENCE" == "high" ]] || [[ "\$CONFIDENCE" == "medium" ]] || [[ \$CURRENT_READS -ge \$MAX_READS ]]; then
            break
        fi
        
        echo "Low confidence (\$CONFIDENCE)"
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
    references = Channel.fromPath("${params.reference_dir}/*.fasta").collect()
    blast_db = CREATE_BLAST_DB(references)
    
    SEROTYPE_FROM_READS(reads_ch, blast_db.blast_db)
    
    emit:
    serotype_info = SEROTYPE_FROM_READS.out.serotype_info
    blast_results = SEROTYPE_FROM_READS.out.blast_results
    blast_db      = blast_db.blast_db
}