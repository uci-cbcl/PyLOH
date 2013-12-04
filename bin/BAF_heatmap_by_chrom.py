#!/usr/bin/env python

import numpy as np
import scipy
import sys
import os
from matplotlib import pyplot
from joint_snv_mix.file_formats.jcnt import JointCountsReader

counts_min = 10
counts_max = 95 
chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
              'chrX', 'chrY']

def main():
    infile = sys.argv[1]
    outdir = sys.argv[2]
    
    os.mkdir(outdir)
    
    reader = JointCountsReader(infile)
    
    chrom_set = reader.get_table_list()
    
    Bins = np.array(range(0, 100 + 1))/100.0
    
    for chrom in chrom_list:
        if chrom not in chrom_set:
            continue
        
        print 'plot %s...' % (chrom)
        
        counts = reader.get_counts(chrom)*1.0
        
        BAF_N = counts[:, 1]/(counts[:, 0] + counts[:, 1])
        BAF_T = counts[:, 3]/(counts[:, 2] + counts[:, 3])
    
        H, _, _ = np.histogram2d(BAF_N, BAF_T, bins=(Bins, Bins))
        
        H_sub = H[counts_min:counts_max, counts_min:counts_max]
        
        color_max = H_sub.max()
    
        pyplot.figure(figsize=(8,8))
        pyplot.xlim((0, 100))
        pyplot.ylim((0, 100))
        pyplot.xticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
        pyplot.yticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
        pyplot.xlabel('Tumour genome B allele frequency (BAF)')
        pyplot.ylabel('Normal genome B allele frequency (BAF)')
        pyplot.imshow(H, vmin = 0, vmax = color_max)#, interpolation='bicubic')
        pyplot.colorbar(ticks=[0, color_max], orientation='vertical')
        pyplot.savefig('./' + outdir + '/' + chrom)
    
    reader.close()

if __name__ == '__main__':
    main()
