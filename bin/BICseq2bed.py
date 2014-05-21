#!/usr/bin/env python

#=======================================================================================================================
# Created on 2013-09-28
# @author: Yi Li
#
# PyLOH
# Copyright (c) 2013 Yi Li <yil8@uci.edu>
#
# This code is free software; you can redistribute it and/or modify it
# under the terms of GNU GPL v2.0 (see the file LICENSE included with the distribution).
#=======================================================================================================================
import argparse

CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
              '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

def main():
    parser = argparse.ArgumentParser(description='Converting BICseq segments file to BED file')

    parser.add_argument('inBICseq', help='''Input BICseq file of segments.''')
    parser.add_argument('outBED', help='''Output bed file of segments.''')
    parser.add_argument('--seg_length', default=1000000, type=int,
                            help='''Minimum length (bp) required for each segment. Default is 1e6.''')
    
    args = parser.parse_args()
    
    infile = open(args.inBICseq)
    outfile = open(args.outBED, 'w')
    
    for line in infile:
        if line[0:5] == 'chrom':
            continue
        
        chrom, start, end = line.strip('\n').split('\t')[0:3]
        
        seg_length = int(end) - int(start) + 1
        
        if seg_length < args.seg_length or chrom not in CHROM_LIST:
            continue
        
        outfile.write('\t'.join([chrom, start, end]) + '\n')
        
    
    infile.close()
    outfile.close()

if __name__ == '__main__':
    main()
