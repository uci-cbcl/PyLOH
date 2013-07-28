'''
Created on 2013-07-27

@author: Yi Li
'''
import numpy as np

from pyloh import constants

def BEDParser(bed_file_name):
    inbed = open(bed_file_name)
    
    chroms = []
    starts = []
    ends = []
    
    for line in inbed:
        if line[0:3] != 'chr':
            continue
        
        chrom, start, end = line.split('\t')[0:3]
        
        chroms.append(chrom)
        starts.append(int(start))
        ends.append(int(end))
            
    inbed.close()
    
    return (chroms, starts, ends)
    
def normal_heterozygous_filter(counts):
    BAF_N_MAX = constants.BAF_N_MAX
    BAF_N_MIN = constants.BAF_N_MIN
    
    I = counts.shape[0]
    idx_keep = []
    
    for i in xrange(0, I):
        a_N = counts[i, 0]
        b_N = counts[i, 1]
        d_N = a_N + b_N
        BAF_N = b_N/d_N
        
        if BAF_N >= BAF_N_MIN and BAF_N <= BAF_N_MAX:
            idx_keep.append(i)
            
    counts = counts[idx_keep]
    
    return counts