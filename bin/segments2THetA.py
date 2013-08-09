#!/usr/bin/env python

import numpy as np
import sys

def main():
    infile = open(sys.argv[1])
    outfile = open(sys.argv[2], 'w')
    
    for line in infile:
        if line[0] == '#':
            outfile.write('#ID' + '\n')
            continue
        
        seg_name, chrom, start, end, normal_reads_num, tumor_reads_num = line.strip('\n').split('\t')[0:6]
        chrom = chrom.strip('chr')
        outfile.write('\t'.join([seg_name, chrom, start, end, tumor_reads_num, normal_reads_num]) + '\n')
        
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main()