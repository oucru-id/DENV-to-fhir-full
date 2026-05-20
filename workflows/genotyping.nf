#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GENOTYPE_NEXTCLADE {
    publishDir "${params.results_dir}/genotyping", mode: 'copy'

    input:
    tuple val(sample_id), path(consensus), path(serotype_json)

    output:
    tuple val(sample_id), path("${sample_id}.genotype_lineage.json"), emit: genotype_info
    path "${sample_id}.nextclade.csv", emit: csv_report
    tuple val(sample_id), path("${sample_id}.nextclade.csv"), emit: raw_nextclade

    script:
    """
    python3 -c "
import json
try:
    with open('${serotype_json}') as f:
        data = json.load(f)
        s = data.get('serotype', '').upper()
        
        if 'DENV-1' in s or 'DENV1' in s:
            print('community/v-gen-lab/dengue/denv1')
        elif 'DENV-2' in s or 'DENV2' in s:
            print('community/v-gen-lab/dengue/denv2')
        elif 'DENV-3' in s or 'DENV3' in s:
            print('community/v-gen-lab/dengue/denv3')
        elif 'DENV-4' in s or 'DENV4' in s:
            print('community/v-gen-lab/dengue/denv4')
        elif 'DENV' in s: 
            print('nextstrain/dengue/all')
        else: 
            print('unknown')
except:
    print('unknown')
" > dataset_name.txt

    DATASET=\$(cat dataset_name.txt)

    if [ "\$DATASET" == "unknown" ]; then
        echo "Unknown serotype, skipping Nextclade"
        echo '{"sample_id": "${sample_id}", "genotype": "unknown", "clade": "unknown", "mutations": []}' > ${sample_id}.genotype_lineage.json
        touch ${sample_id}.nextclade.csv
    else
        nextclade dataset get --name \$DATASET --output-dir nextclade_db

        nextclade run \\
            --input-dataset nextclade_db \\
            --output-json ${sample_id}.nextclade.json \\
            --output-csv ${sample_id}.nextclade.csv \\
            ${consensus}

        python3 -c "
import json
import csv
import re
import sys

try:
    with open('${serotype_json}') as f:
        sero_data = json.load(f)
        real_serotype = sero_data.get('serotype', 'Unknown')

    csv_file = '${sample_id}.nextclade.csv'
    results = {}
    
    with open(csv_file, newline='') as f:
        reader = csv.DictReader(f, delimiter=';')
        for row in reader:
            results = row
            break
    
    if results:
        clade = results.get('clade', 'unknown') # e.g., '3III_B.3.2' or '1I.A'
        aa_subs_str = results.get('aaSubstitutions', '')
        qc_score = results.get('qc.overallScore', 0)
        coverage = results.get('coverage', 0)

        genotype = 'Unclassified'
        major_lineage = 'unknown'
        minor_lineage = 'unknown'
        
        match = re.match(r'^(\\d)([IVX]+)(?:[._]([A-Z0-9]+))?(?:[._](.+))?\$', clade)
        
        if match:
            parsed_geno = match.group(2) # 'I', 'II', 'III'
            parsed_maj  = match.group(3) # 'A', 'B'
            parsed_min  = match.group(4) # '1.2', '3.2'
            
            genotype = f'Genotype {parsed_geno}'
            if parsed_maj: major_lineage = f'Lineage {parsed_maj}'
            if parsed_min: minor_lineage = f'Sub-lineage {parsed_min}'
        
        elif clade.startswith('DENV'):
             genotype = f'{clade} (Basal/Unclassified)'
        
        mutations = []
        if aa_subs_str:
            mutations = [m.strip() for m in aa_subs_str.split(',') if m.strip()]

        output = {
            'sample_id': '${sample_id}',
            'serotype': real_serotype, 
            'genotype': genotype,
            'major_lineage': major_lineage, 
            'minor_lineage': minor_lineage,
            'clade': clade,           
            'confidence': 'high',
            'method': 'nextclade',
            'mutations': mutations,
            'analysis': {
                'qc_score': float(qc_score) if qc_score else 0,
                'coverage': float(coverage) if coverage else 0
            }
        }
    else:
        output = {'sample_id': '${sample_id}', 'genotype': 'unknown', 'error': 'No results in CSV'}

except Exception as e:
    output = {
        'sample_id': '${sample_id}',
        'genotype': 'error',
        'error': str(e)
    }

with open('${sample_id}.genotype_lineage.json', 'w') as f:
    json.dump(output, f, indent=2)
"
    fi
    """
}

workflow GENOTYPING {
    take:
    consensus_ch
    serotype_ch

    main:
    input_ch = consensus_ch.join(serotype_ch)
    
    GENOTYPE_NEXTCLADE(input_ch)
    
    emit:
    genotype_info = GENOTYPE_NEXTCLADE.out.genotype_info  
    csv_report    = GENOTYPE_NEXTCLADE.out.csv_report
    raw_nextclade = GENOTYPE_NEXTCLADE.out.raw_nextclade
}
