#!/usr/bin/env python

import numpy as np
import sys
import pysam

CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

def main():
    innormal = sys.argv[1]
    intumor = sys.argv[2]
    outfile = sys.argv[3]
    
    normal_bam = pysam.Samfile(innormal,'rb')
    tumor_bam = pysam.Samfile(intumor,'rb')
    outfile = open(outfile, 'w')
    
    references = normal_bam.references
    lengths = normal_bam.lengths
    ref_num = len(references)
    chrom_len = {}
    
    print '\t'.join(['#ID', 'chrm', 'start', 'end', 'tumorCount', 'normalCount'])
    outfile.write('\t'.join(['#ID', 'chrm', 'start', 'end', 'tumorCount', 'normalCount']) + '\n')
    
    for i in range(0, ref_num):
        if references[i] not in CHROM_LIST:
            continue
                
        chrom = references[i].strip('chr')
        start = 0
        end = lengths[i] - 1
        ID = 'start_%s_%s:end_%s_%s' % (chrom, start, chrom, end)
        d_N = normal_bam.count(references[i])
        d_T = tumor_bam.count(references[i])
        
        print '\t'.join([ID, chrom, str(start), str(end), str(d_T), str(d_N)])
        outfile.write('\t'.join([ID, chrom, str(start), str(end), str(d_T), str(d_N)]) + '\n')
        
    normal_bam.close()
    tumor_bam.close()
    outfile.close()

if __name__ == '__main__':
    main()