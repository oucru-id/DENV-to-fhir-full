#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process REMOVE_HOST_HOSTILE {
    publishDir "${params.results_dir}/host_removal", mode: 'copy', pattern: "*.hostile_log.txt"
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    tuple val(sample_id), path("${sample_id}_clean*.fastq.gz"), emit: host_removed
    path("${sample_id}.hostile_log.txt"), emit: log
    
    script:
    if (reads instanceof List) {
        """
        echo "Processing paired-end reads for ${sample_id}"
        
        hostile clean \\
            --fastq1 ${reads[0]} \\
            --fastq2 ${reads[1]} \\
            --aligner bowtie2 \\
            --threads ${task.cpus} \\
            --output . \\
            2>&1 | tee ${sample_id}.hostile_log.txt
        
        for file in *.clean*.fastq.gz; do
            if [[ \$file == *".clean_1.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_1.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_1.fastq.gz"
            elif [[ \$file == *".clean_2.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_2.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_2.fastq.gz"
            elif [[ \$file == *"_1.clean.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_1.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_1.fastq.gz"
            elif [[ \$file == *"_2.clean.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_2.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_2.fastq.gz"
            elif [[ \$file == *".clean.1.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_1.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_1.fastq.gz"
            elif [[ \$file == *".clean.2.fastq.gz" ]]; then
                mv "\$file" "${sample_id}_clean_2.fastq.gz"
                echo "Renamed \$file to ${sample_id}_clean_2.fastq.gz"
            fi
        done
        
        if [[ ! -f "${sample_id}_clean_1.fastq.gz" ]] || [[ ! -f "${sample_id}_clean_2.fastq.gz" ]]; then
            exit 1
        fi

        """
    } else {
        """
        echo "Processing single-end reads for ${sample_id}"
        
        hostile clean \\
            --fastq1 ${reads} \\
            --aligner minimap2 \\
            --threads ${task.cpus} \\
            --output . \\
            2>&1 | tee ${sample_id}.hostile_log.txt
        
        found=false
        for file in *.clean*.fastq.gz; do
            mv "\$file" "${sample_id}_clean.fastq.gz"
            echo "Renamed \$file to ${sample_id}_clean.fastq.gz"
            found=true
            break
        done
        
        if [[ ! -f "${sample_id}_clean.fastq.gz" ]]; then
            ls -la *.fastq.gz
            exit 1
        fi
        
        echo "Final output file:"
        ls -la ${sample_id}_clean*.fastq.gz
        """
    }
}

workflow HOST_REMOVAL {
    take:
    reads_ch
    sequencing_type  
    
    main:
    host_removed = REMOVE_HOST_HOSTILE(reads_ch)
    
    emit:
    host_removed = host_removed.host_removed

}
