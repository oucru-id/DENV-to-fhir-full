#!/usr/bin/env python3

import argparse
import json
import pandas as pd
import re
import sys
import traceback

def classify_serotype(blast_file, sample_id, total_input_reads):
    try:
        try:
            df = pd.read_csv(blast_file, sep='\t', header=None,
                            names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                   'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                                   'evalue', 'bitscore', 'stitle'])
        except pd.errors.EmptyDataError:
            df = pd.DataFrame()

        if df.empty:
            return {
                'sample_id': sample_id,
                'serotype': 'unknown',
                'confidence': 'very_low',
                'strain_type': 'unknown',
                'top_match': None,
                'identity': 0,
                'total_input_reads': total_input_reads,
                'mapped_reads': 0,
                'note': 'No BLAST hits found'
            }
        
        df = df[df['pident'] > 75]
        
        if df.empty:
             return {
                'sample_id': sample_id,
                'serotype': 'unknown',
                'confidence': 'very_low',
                'total_input_reads': total_input_reads,
                'mapped_reads': 0,
                'note': 'No hits passed quality filter (>75% identity)'
            }

        df = df.sort_values(['qseqid', 'bitscore'], ascending=[True, False])
        best_hits = df.drop_duplicates('qseqid', keep='first')
        mapped_reads = len(best_hits)
        serotype_counts = {}
        serotype_stats = {} 

        for _, hit in best_hits.iterrows():
            title = str(hit['stitle']).upper()
            serotype_match = re.search(r'DENGUE\s+(?:VIRUS)?\s*(?:TYPE)?\s*([1-4])', title)
            
            if serotype_match:
                serotype_num = serotype_match.group(1)
                serotype_key = f"DENV-{serotype_num}"
                
                if serotype_key not in serotype_counts:
                    serotype_counts[serotype_key] = 0
                    serotype_stats[serotype_key] = {'pident_sum': 0, 'bitscore_sum': 0, 'count': 0}
                
                serotype_counts[serotype_key] += 1
                serotype_stats[serotype_key]['pident_sum'] += hit['pident']
                serotype_stats[serotype_key]['bitscore_sum'] += hit['bitscore']
                serotype_stats[serotype_key]['count'] += 1

        if not serotype_counts:
            return {
                'sample_id': sample_id,
                'serotype': 'unknown',
                'confidence': 'very_low',
                'total_input_reads': total_input_reads,
                'mapped_reads': mapped_reads,
                'note': 'Hits found but could not parse'
            }

        sorted_serotypes = sorted(serotype_counts.items(), key=lambda x: x[1], reverse=True)
        winner_serotype, winner_count = sorted_serotypes[0]
        total_classified_reads = sum(serotype_counts.values())
        winner_proportion = winner_count / total_classified_reads
        confidence = 'low'
        confidence_penalty = False
        consistency = 'consistent'
        
        if len(sorted_serotypes) > 1:
            second_serotype, second_count = sorted_serotypes[1]
            second_proportion = second_count / total_classified_reads
            
            if second_proportion > 0.10:
                confidence_penalty = True
                consistency = 'mixed'

        avg_pident = serotype_stats[winner_serotype]['pident_sum'] / serotype_stats[winner_serotype]['count']
        
        if winner_count > 50 and avg_pident > 95 and winner_proportion > 0.9:
            confidence = 'high'
        elif winner_count > 10 and avg_pident > 90 and winner_proportion > 0.7:
            confidence = 'medium'
        elif winner_count < 5:
            confidence = 'very_low'
        
        if confidence_penalty:
            if confidence == 'high': confidence = 'medium'
            elif confidence == 'medium': confidence = 'low'

        winner_hits = best_hits[best_hits['stitle'].str.contains(f"type {winner_serotype.split('-')[1]}", case=False, regex=False)]
        if winner_hits.empty:
             winner_hits = best_hits
             
        top_match = winner_hits.sort_values('bitscore', ascending=False).iloc[0]
        ref_id = top_match['sseqid']
        
        strain_type = 'human'
        if any(x in str(ref_id).upper() for x in ['SYLVATIC']):
            strain_type = 'sylvatic'

        return {
            'sample_id': sample_id,
            'serotype': winner_serotype,
            'strain_type': strain_type,
            'confidence': confidence,
            'top_match': ref_id,
            'identity': float(avg_pident),
            'support_reads': int(winner_count),
            'mapped_reads': int(mapped_reads),
            'total_input_reads': int(total_input_reads),
            'serotype_consistency': consistency,
            'analysis_details': {
                'serotype_counts': serotype_counts,
                'serotype_proportions': {k: round(v/total_classified_reads, 4) for k,v in serotype_counts.items()},
                'winner_proportion': round(winner_proportion, 4)
            }
        }
        
    except Exception as e:
        return {
            'sample_id': sample_id,
            'serotype': 'unknown',
            'strain_type': 'unknown',
            'confidence': 'very_low',
            'error': str(e),
            'error_trace': traceback.format_exc()
        }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--blast_results', required=True)
    parser.add_argument('--sample_id', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--total_reads', required=False, default=0, type=int, help="Total reads in input FASTQ")
    args = parser.parse_args()
    
    result = classify_serotype(args.blast_results, args.sample_id, args.total_reads)
    
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"Serotype classification: {result['serotype']} (confidence: {result.get('confidence', 'unknown')})")
    
    if 'mapped_reads' in result:
        print(f"  Mapped reads: {result['mapped_reads']} / {result.get('total_input_reads', '?')}")
        
    if 'analysis_details' in result:
        details = result['analysis_details']
        if 'serotype_counts' in details:
            print(f"  Serotype counts: {details['serotype_counts']}")

if __name__ == '__main__':
    main()