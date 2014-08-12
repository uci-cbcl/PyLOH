#!/usr/bin/env python

#=======================================================================================================================
# Created on 2014-08-11
# @author: Yi Li
#
# PyLOH
# Copyright (c) 2014 Yi Li <yil8@uci.edu>
#
# This code is free software; you can redistribute it and/or modify it
# under the terms of GNU GPL v2.0 (see the file LICENSE included with the distribution).
#=======================================================================================================================
import argparse

import sys

import numpy as np
import pysam


CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
              '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
              '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']


def main():
    parser = argparse.ArgumentParser(description='Converting paired normal and tumor BAM files to input for DNAcopy')

    parser.add_argument('normal_bam', help='''Input BAM file for normal sample.''')
    parser.add_argument('tumor_bam', help='''Input BAM file for tumor sample.''')
    parser.add_argument('exons_bed', help='''Input BED file for all exon regions.''')
    parser.add_argument('DNAcopy_bed', help='''Output BED file for DNAcopy.''')
    parser.add_argument('--min_depth', default=100, type=int,
                            help='''Minimum reads detph required for each exon region in both
                            normal and tumor samples. Default is 100.''')
    
    args = parser.parse_args()
    
    
    
    normal_bam = pysam.Samfile(args.normal_bam, 'rb')
    tumor_bam = pysam.Samfile(args.tumor_bam, 'rb')
    exons_bed = open(args.exons_bed)
    DNAcopy_bed = open(args.DNAcopy_bed, 'w')
    depth_thred = args.min_depth
    
    i = 0
    for line in exons_bed:
        i+= 1
        if i % 1000 == 0:
            print '%s exons processed...' % (i)
            sys.stdout.flush()

        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        
        normal_reads = normal_bam.count(chrom, start, end)
        tumor_reads = tumor_bam.count(chrom, start, end)
        
        if normal_reads < args.min_depth or tumor_reads < args.min_depth or chrom not in CHROM_LIST:
            continue
        
        log2_ratio = np.log2(tumor_reads*1.0/normal_reads)
        pos = (end + start)/2
        
        fields.extend(map(str, [normal_reads, tumor_reads, log2_ratio, pos]))
        
        DNAcopy_bed.write('\t'.join(fields) + '\n')
    
    normal_bam.close()
    tumor_bam.close()
    exons_bed.close()
    DNAcopy_bed.close()


if __name__ == '__main__':
    main()
