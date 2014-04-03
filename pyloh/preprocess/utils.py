'''
Created on 2013-07-27

@author: Yi Li
'''
import sys
import numpy as np

from pyloh import constants

def BEDParser(bed_file_name):
    inbed = open(bed_file_name)
    
    chroms = []
    starts = []
    ends = []
    
    for line in inbed:
        fields = line.split('\t')
        chrom_name = fields[0]
        chrom_ID = chrom_name_to_ID(chrom_name)
        
        if chrom_ID == -1:
            continue
        
        chrom_name, start, end = fields[0:3]
        
        chroms.append(chrom_name)
        starts.append(int(start))
        ends.append(int(end))
            
    inbed.close()
    
    return (chroms, starts, ends)
    
def chrom_ID_to_name(ID, format):
    if format == 'UCSC':
        chrom_name = 'chr' + str(ID)
    elif format == 'ENSEMBL':
        chrom_name = str(ID)
    else:
        print 'Error: %s not supported' % (format)
        sys.exit(1)

    return chrom_name

def chrom_name_to_ID(chrom_name):
    ID = -1
    
    try:
        ID = int(chrom_name.strip('chr'))
    except:
        pass
        
    return ID

def get_chrom_format(chroms):
    format = 'NONE'
    
    for chrom in chroms:
        if chrom[0:3] == 'chr':
            format = 'UCSC'
            break
        else:
            try:
                ID = int(chrom)
                format = 'ENSEMBL'
                break
            except:
                pass
    
    if format == 'NONE':
        print 'Error: %s not supported' % (chrom)
        sys.exit(-1)
    else:
        return format

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

def get_BAF_counts(counts):
    BAF_bins = constants.BAF_BINS
    
    a_N = counts[:, 0]*1.0
    b_N = counts[:, 1]*1.0
    a_T = counts[:, 2]*1.0
    b_T = counts[:, 3]*1.0
    
    BAF_N = b_N/(a_N + b_N)
    BAF_T = b_T/(a_T + b_T)
    
    BAF_counts, _, _ = np.histogram2d(BAF_N, BAF_T, bins=(BAF_bins, BAF_bins))
    
    return BAF_counts

def tumor_LOH_test(counts, WES_flag):
    BAF_T_MAX = constants.BAF_T_MAX
    BAF_T_MIN = constants.BAF_T_MIN
    LOH_FREC_MAX = constants.LOH_FREC_MAX
    LOH_FREC_MIN = constants.LOH_FREC_MIN
    
    I = counts.shape[0]
    
    if WES_flag == True:
        sites_num_min = constants.SITES_NUM_MIN_WES
    else:
        sites_num_min = constants.SITES_NUM_MIN_WGS
    
    if I < sites_num_min:
        LOH_frec = -1
        LOH_status = 'NONE'
        
        return (LOH_frec, LOH_status)
    
    LOH_num = 0.0
    
    for i in xrange(0, I):
        a_T = counts[i, 2]*1.0
        b_T = counts[i, 3]*1.0
        d_T = a_T + b_T
        BAF_T = b_T/d_T
        
        if BAF_T < BAF_T_MIN or BAF_T > BAF_T_MAX:
            LOH_num = LOH_num + 1
    
    LOH_frec = LOH_num/I
    
    if LOH_frec < LOH_FREC_MIN:
        LOH_status = 'FALSE'
    elif LOH_frec >= LOH_FREC_MIN and LOH_frec < LOH_FREC_MAX:
        LOH_status = 'UNCERTAIN'
    elif LOH_frec >= LOH_FREC_MAX:
        LOH_status = 'TRUE'
    else:
        LOH_status = 'ERROR'
        
    return (LOH_frec, LOH_status)
    
def remove_outliers(X):
    idx_keep = []
    
    n = X.shape[0]
    
    for i in range(0, n):
        if np.abs(X[i] - X.mean()) <= X.std():
            idx_keep.append(i)
            
    X = X[idx_keep]
    
    return X
    