#!/usr/bin/env python3

import argparse
import json
import os
from datetime import datetime
from Bio import SeqIO

def load_json_file(file_path):
    try:
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            return {}
        with open(file_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {file_path}: {e}")
        return {}

def get_sequence_stats(fasta_file):
    try:
        if not os.path.exists(fasta_file):
            return {"length": 0, "n_count": 0, "gc_content": 0}
        
        with open(fasta_file, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                seq = str(record.seq).upper()
                length = len(seq)
                n_count = seq.count('N')
                gc_count = seq.count('G') + seq.count('C')
                gc_content = (gc_count / length * 100) if length > 0 else 0
                
                return {
                    "length": length,
                    "n_count": n_count,
                    "n_percentage": (n_count / length * 100) if length > 0 else 0,
                    "gc_content": round(gc_content, 2),
                    "coverage": round((length - n_count) / length * 100, 2) if length > 0 else 0
                }
        
        return {"length": 0, "n_count": 0, "gc_content": 0, "coverage": 0}
    except Exception as e:
        return {"length": 0, "n_count": 0, "gc_content": 0, "coverage": 0}

def generate_dengue_report(sample_id, serotype_data, genotype_data, consensus_file, output_file):
    seq_stats = get_sequence_stats(consensus_file)
    
    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write(f"DENGUE VIRUS GENOMIC ANALYSIS REPORT\n")
        f.write(f"Sample ID: {sample_id}\n")
        f.write("="*80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        f.write("SEQUENCE QUALITY SUMMARY:\n")
        f.write("-" * 50 + "\n")
        f.write(f"Consensus length: {seq_stats['length']:,} bp\n")
        f.write(f"Genome coverage: {seq_stats['coverage']:.1f}%\n")
        f.write(f"Ambiguous bases (N): {seq_stats['n_count']:,} ({seq_stats['n_percentage']:.1f}%)\n")
        f.write(f"GC content: {seq_stats['gc_content']:.1f}%\n\n")
        
        f.write("SEROTYPE CLASSIFICATION:\n")
        f.write("-" * 50 + "\n")
        serotype = serotype_data.get('serotype', 'Unknown')
        
        f.write(f"Serotype: {serotype}\n")
        f.write("\n")
        
        f.write("GENOTYPE CLASSIFICATION:\n")
        f.write("-" * 50 + "\n")
        
        if genotype_data:
            genotype = genotype_data.get('genotype', 'Unknown')
            major_lineage = genotype_data.get('major_lineage', 'Unknown')
            minor_lineage = genotype_data.get('minor_lineage', 'Unknown')
            clade = genotype_data.get('clade', 'Unknown')
            
            f.write(f"Genotype: {genotype}\n")
            f.write(f"Major lineage: {major_lineage}\n")
            f.write(f"Minor lineage: {minor_lineage}\n")
            f.write(f"Clade Code: {clade}\n")
            
            mutations = genotype_data.get('mutations', [])
            if mutations:
                f.write(f"\nDetected Amino Acid Mutations ({len(mutations)}):\n")
                from collections import defaultdict
                genes = defaultdict(list)
                for m in mutations:
                    if ':' in m:
                        g, change = m.split(':')
                        genes[g].append(change)
                
                for gene, changes in sorted(genes.items()):
                    f.write(f"  {gene}: {', '.join(changes)}\n")
            
        else:
            f.write("No genotyping data available\n")
        f.write("\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("END OF DENGUE ANALYSIS REPORT\n")
        f.write("="*80 + "\n")

def main():
    parser = argparse.ArgumentParser(description='Generate dengue virus analysis report')
    parser.add_argument('--sample_id', required=True, help='Sample identifier')
    parser.add_argument('--serotype_json', required=True, help='Serotype JSON file')
    parser.add_argument('--genotype_json', help='Genotype JSON file (optional)')
    parser.add_argument('--consensus', required=True, help='Consensus sequence FASTA file')
    parser.add_argument('--output', required=True, help='Output report file')
    
    args = parser.parse_args()
    
    serotype_data = load_json_file(args.serotype_json)
    genotype_data = load_json_file(args.genotype_json) if args.genotype_json else {}
    
    generate_dengue_report(
        args.sample_id,
        serotype_data,
        genotype_data,
        args.consensus,
        args.output
    )
    
if __name__ == '__main__':
    main()