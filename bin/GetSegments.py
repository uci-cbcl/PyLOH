#!/usr/bin/env python


centromere_starts = {}
centromere_ends = {}
centromere_starts['chr1'] = 121100000
centromere_ends['chr1'] = 128000000
chrom_start = 0

chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

chrom_lens = [247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
              158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
              114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
              63811651, 62435964, 46944323, 49691432]

import sys
import numpy as np

def main():
    outfile = open(sys.argv[1], 'w')
    chrom = sys.argv[2]
    p_seg_num = int(sys.argv[3])
    q_seg_num = int(sys.argv[4])
    
    chrom_idx = chrom_list.index(chrom)
    p_start = chrom_start
    p_end = centromere_starts[chrom]
    q_start = centromere_ends[chrom]
    q_end = chrom_lens[chrom_idx]

    p_alpha = np.random.random_integers(1, p_seg_num, p_seg_num)
    q_alpha = np.random.random_integers(1, q_seg_num, q_seg_num)
    
    p_percent = p_alpha*1.0/p_alpha.sum()
    q_percent = q_alpha*1.0/q_alpha.sum()
    
    p_seg_lens = np.floor((p_end - p_start)*p_percent)
    q_seg_lens = np.floor((q_end - q_start)*q_percent)
    
    outfile.write('\t'.join(['#chrom', 'start', 'end', 'paternal_copy', 'maternal_copy']) + '\n')
    
    base = p_start
    for i in range(0, p_seg_num):
        start = int(base)
        end = int(base + p_seg_lens[i])
        outfile.write('\t'.join([chrom, str(start), str(end)]) + '\n')
        base = base + p_seg_lens[i]
        
    base = q_start
    for i in range(0, q_seg_num):
        start = int(base)
        end = int(base + q_seg_lens[i])
        outfile.write('\t'.join([chrom, str(start), str(end)]) + '\n')
        base = base + q_seg_lens[i]
    
    outfile.close()
    
if __name__ == '__main__':
    main()