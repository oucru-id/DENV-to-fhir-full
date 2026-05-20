#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process fastqc {
    
    input:
    tuple val(sample_id), path(reads)  

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: qc_report

    script:
    """
    fastqc ${reads.join(' ')} --outdir .
    """
}

process select_reference {
    
    input:
    tuple val(sample_id), path(reads), path(serotype_json)
    
    output:
    tuple val(sample_id), path(reads), path("selected_reference.fasta"), emit: reads_with_ref
    
    script:
    """
    #!/usr/bin/env python3
    import json
    import shutil
    import os
    import glob
    
    try:
        with open('${serotype_json}', 'r') as f:
            data = json.load(f)
        
        if 'classification' in data:
            serotype = data['classification'].get('serotype', 'unknown')
        else:
            serotype = data.get('serotype', 'unknown')
            
        if serotype and serotype != 'unknown':
            serotype = serotype.replace('-', '').upper()
        
        if serotype == 'unknown' and 'top_match' in data:
            top_match = data['top_match']
            print(f"Serotype unknown, inferring from top_match: {top_match}")
            
            if 'NC_001477' in top_match or 'DENV-1' in top_match or 'DENV1' in top_match:
                serotype = 'DENV1'
            elif 'NC_001474' in top_match or 'DENV-2' in top_match or 'DENV2' in top_match:
                serotype = 'DENV2'
            elif 'NC_001475' in top_match or 'DENV-3' in top_match or 'DENV3' in top_match:
                serotype = 'DENV3'
            elif 'NC_002640' in top_match or 'DENV-4' in top_match or 'DENV4' in top_match:
                serotype = 'DENV4'
            elif 'EF105380' in top_match:
                serotype = 'SYLVATIC1'
            elif 'JF262780' in top_match:
                serotype = 'SYLVATIC2'
            
            print(f"serotype: {serotype}")
    except Exception as e:
        print(f"Error parsing JSON: {e}")
        serotype = 'DENV1'
    
    ref_map = {
        'DENV1': 'NC_0014771.fasta',    
        'DENV2': 'NC_0014742.fasta',    
        'DENV3': 'NC_0014752.fasta',    
        'DENV4': 'NC_0026401.fasta',   
        'SYLVATIC1': 'EF105380.fasta',  
        'SYLVATIC2': 'JF262780.fasta'   
    }
    
    ref_file = ref_map.get(serotype, 'NC_0014771.fasta')
    print(f"Selected serotype: {serotype}, mapping to file: {ref_file}")

    ref_path = os.path.join('${params.reference_dir}', ref_file)
    
    if not os.path.exists(ref_path):
        all_fastas = glob.glob(os.path.join('${params.reference_dir}', '*.fasta'))
        if all_fastas:
            ref_path = all_fastas[0]
            print(f"WARNING: Requested reference not found, falling back to: {ref_path}")
        else:
            raise FileNotFoundError(f"No FASTA reference files found in ${params.reference_dir}")
    
    print(f"Using reference: {os.path.basename(ref_path)}")
    
    shutil.copy(ref_path, 'selected_reference.fasta')
    """
}

process bwa_mem2 {

    input:
    tuple val(sample_id), path(reads), path(reference)

    output:
    tuple val(sample_id), path("aligned.bam"), path(reference), emit: aligned

    script:
    """
    if [ ! -f ${reference}.bwt.2bit.64 ]; then
        bwa-mem2 index ${reference}
    fi

    bwa-mem2 mem -t ${task.cpus} \
        -R "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:lib1" \
        ${reference} \
        ${reads[0]} ${reads[1]} 2> bwa.log | \
    samtools sort -@ ${task.cpus} -m 2G -o aligned.bam -

    if [ ! -s aligned.bam ]; then
        echo "ERROR: BAM file is empty"
        cat bwa.log
        exit 1
    fi

    samtools index aligned.bam
    
    samtools flagstat aligned.bam > flagstat.txt
    cat flagstat.txt
    """
}

process generate_consensus {

    publishDir "${params.results_dir}/consensus", mode: 'copy'

    input:
    tuple val(sample_id), path("aligned.bam"), path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fasta"),                        emit: consensus
    tuple val(sample_id), path("${sample_id}.vcf.gz"), path("${sample_id}.vcf.gz.tbi"), emit: vcf
    tuple val(sample_id), path("${sample_id}.coverage.txt"),                          emit: coverage_stats

    script:
    """
    if [ ! -f aligned.bam.bai ]; then
        samtools index aligned.bam
    fi

    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi

    bcftools mpileup -f ${reference} aligned.bam | \
    bcftools call -mv -Oz -o ${sample_id}.vcf.gz

    tabix -p vcf ${sample_id}.vcf.gz

    bcftools consensus -f ${reference} ${sample_id}.vcf.gz > temp_consensus.fasta

    sed "s/^>.*/>${sample_id}_consensus/" temp_consensus.fasta > ${sample_id}_consensus.fasta

    samtools coverage aligned.bam > ${sample_id}.coverage.txt
    """
}

workflow ILLUMINA {
    take:
    reads_with_serotype

    main:
    reads_only = reads_with_serotype.map { id, reads, json -> tuple(id, reads) }
    qc_clean = fastqc(reads_only)
    
    reads_with_ref = select_reference(reads_with_serotype)
    
    aligned = bwa_mem2(reads_with_ref.reads_with_ref)
    consensus_out = generate_consensus(aligned.aligned)

    emit:
    qc_report = qc_clean.qc_report
    consensus = consensus_out.consensus
    vcf       = consensus_out.vcf
    coverage  = consensus_out.coverage_stats
}
