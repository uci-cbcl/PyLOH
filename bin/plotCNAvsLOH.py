#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
import scipy.stats as st
from matplotlib import pyplot as plt

from pyloh.preprocess.data import Data, Segments

MU_0 = 0.5
P_VALUE = 0.05



def ab2ld(counts):
    l = np.minimum(counts[:, 2], counts[:, 3])
    d = counts[:, 2] + counts[:, 3]
    
    return l, d

def getLOHfrec(counts):
    l, d = ab2ld(counts)
    
    p_T = st.binom.cdf(l, d, MU_0)
    LOH_num = np.where(p_T < P_VALUE)[0].shape[0]
    sites_num = p_T.shape[0]
    
    return LOH_num*1.0/sites_num
    

def main():
    filename_base = sys.argv[1]
    
    data = Data()
    data.read_data(filename_base)
    
    paired_counts = data.paired_counts
    segments = data.segments
    seg_num = data.seg_num
    
    log2_ratio_lst = []
    LOH_frec_lst = []
    for j in range(0, seg_num):
        log2_ratio = segments[j][8]
        LOH_frec = getLOHfrec(paired_counts[j])
        
        log2_ratio_lst.append(log2_ratio)
        LOH_frec_lst.append(LOH_frec)
        
    plt.plot(log2_ratio_lst, LOH_frec_lst, 'bx', linewidth = 5)
    plt.xlim(-2, 2)
    plt.ylim(0, 1)
    plt.xlabel('copy number log2 ratio')
    plt.ylabel('LOH sites fraction')
    plt.show()


if __name__ == '__main__':
    main()