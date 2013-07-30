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
        a_N = counts[i, 0]*1.0
        b_N = counts[i, 1]*1.0
        d_N = a_N + b_N
        BAF_N = b_N/d_N
        
        if BAF_N >= BAF_N_MIN and BAF_N <= BAF_N_MAX:
            idx_keep.append(i)
            
    counts = counts[idx_keep]
    
    return counts

def tumor_LOH_test(counts):
    BAF_T_MAX = constants.BAF_T_MAX
    BAF_T_MIN = constants.BAF_T_MIN
    LOH_FREC_THRED = constants.LOH_FREC_THRED
    
    I = counts.shape[0]
    
    LOH_num = 0.0
    
    for i in xrange(0, I):
        a_T = counts[i, 2]*1.0
        b_T = counts[i, 3]*1.0
        d_T = a_T + b_T
        BAF_T = b_T/d_T
        
        if BAF_T < BAF_T_MIN or BAF_T > BAF_T_MAX:
            LOH_num = LOH_num + 1
    
    LOH_frec = LOH_num/I
    
    if LOH_frec < LOH_FREC_THRED:
        LOH_flag = False
    else:
        LOH_flag = True
        
    return (LOH_frec, LOH_flag)