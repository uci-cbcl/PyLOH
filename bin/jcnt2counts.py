#!/usr/bin/env python

import numpy as np
import scipy
import sys
from matplotlib import pyplot
from joint_snv_mix.file_formats.jcnt import JointCountsReader

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    
    reader = JointCountsReader(infile)
    outNbaf = open(outfile + '.N.baf.txt', 'w')
    outTbaf = open(outfile + '.T.baf.txt', 'w')
    
    outNbaf.write('\t'.join(['chr', 'position', 'coverage', 'baf']) + '\n')
    outTbaf.write('\t'.join(['chr', 'position', 'coverage', 'baf']) + '\n')
    
    chr_list = reader.get_table_list()
    
    for chrom in chr_list:
        jcnt_table = reader.get_table(chrom)
        I = jcnt_table.shape[0]
        for i in xrange(0, I):
            pos, __, __, __, a_N, b_N, a_T, b_T = jcnt_table[i]
            
            outNbaf.write('\t'.join([chrom, str(pos), str(a_N + b_N), str(b_N)]) + '\n')
            outTbaf.write('\t'.join([chrom, str(pos), str(a_T + b_T), str(b_T)]) + '\n')
            
            if (i + 1) % 100000 == 0:
                print '%s %s sites processed...' % (chrom, i + 1)

    reader.close()
    outNbaf.close()
    outTbaf.close()

if __name__ == '__main__':
    main()
