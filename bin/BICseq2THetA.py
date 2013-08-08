#!/usr/bin/env python

import numpy as np
import sys
import pysam

CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

#len_thred = 50000

def main():
    infile = open(sys.argv[1])
    outfile = open(sys.argv[2], 'w')
    len_thred = int(sys.argv[3])
    
    outfile.write('\t'.join(['#ID', 'chrm', 'start', 'end', 'tumorCount', 'normalCount']) + '\n')
    for line in infile:
        if line[0:5] == 'chrom':
            continue
        
        chrom, start, end, tumorCount, normalCount, __, __ = line.strip('\n').split('\t')
        segment_len = int(end) - int(start) + 1
        
        if segment_len < len_thred or chrom not in CHROM_LIST:
            continue
        
        chrom = chrom.strip('chr')
        
        ID = 'start_%s_%s:end_%s_%s' % (chrom, start, chrom, end)

        print('\t'.join([ID, chrom, start, end, tumorCount, normalCount, str(segment_len)]))
        outfile.write('\t'.join([ID, chrom, start, end, tumorCount, normalCount]) + '\n')
        
        
    infile.close()
    outfile.close()

if __name__ == '__main__':
    main()
