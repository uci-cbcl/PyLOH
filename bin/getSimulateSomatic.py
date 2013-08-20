#!/usr/bin/env python

import sys
from collections import Counter

import pysam
import numpy as np

def reads_to_bases(reads, pos):
    bases = []
    
    for read in reads:
        if pos not in read.positions:
            continue
        
        idx = read.positions.index(pos)
        base = read.seq[idx-1:idx].upper()
        bases.append(base)
    
    return bases


def main():
    innormal = sys.argv[1]
    intumor = sys.argv[2]
    infasta = sys.argv[3]
    invcf = sys.argv[4]
    outcounts = sys.argv[5]
    chrom = sys.argv[6]
    
    normal_bam = pysam.Samfile(innormal, 'rb')
    tumor_bam = pysam.Samfile(intumor, 'rb')
    ref_genome_fasta = pysam.Fastafile(infasta)
    invcf = open(invcf)
    outcounts = open(outcounts, 'w')
    
    normal_pileup = normal_bam.pileup(chrom)
    tumor_pileup = tumor_bam.pileup(chrom)
    
    outcounts.write('\t'.join(['#pos', 'ref', 'alt', 'a_N', 'b_N', 'a_T', 'b_T']) + '\n')
    for line in invcf:
        if line[0] == '#':
            continue
        
        chrom_vcf, pos, __, ref, alt = line.strip('\n').split('\t')[0:5]
        pos = int(pos)
        ref = ref.upper()
        alt = alt.upper()
        
        if chrom_vcf != chrom:
            continue
        
        reads_N = normal_bam.fetch(chrom, pos - 1, pos)
        reads_T = tumor_bam.fetch(chrom, pos - 1, pos)
        
        bases_N = reads_to_bases(reads_N, pos)
        bases_T = reads_to_bases(reads_T, pos)
        
        counter_N = Counter(bases_N)
        counter_T = Counter(bases_T)
        
        a_N = counter_N[ref]
        b_N = counter_N[alt]
        d_N = a_N + b_N
        a_T = counter_T[ref]
        b_T = counter_T[alt]
        d_T = a_T + b_T
        
        if b_T <= 2 or b_N >= 2:
            continue
        
        
        outcounts.write('\t'.join(map(str, [pos, ref, alt, a_N, b_N, a_T, b_T])) + '\n')
        
    
    
    normal_bam.close()
    tumor_bam.close()
    ref_genome_fasta.close()
    invcf.close()
    outcounts.close()
    
if __name__ == '__main__':
    main()