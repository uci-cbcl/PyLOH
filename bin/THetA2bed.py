#!/usr/bin/env python

import sys

def main():
    infile = open(sys.argv[1])
    outfile = open(sys.argv[2], 'w')
    
    for line in infile:
        if line[0] == '#':
            continue
        
        chrom_idx, start, end = line.split('\t')[1:4]
        chrom = 'chr' + chrom_idx
        
        outfile.write('\t'.join([chrom, start, end]) + '\n')
        
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main()