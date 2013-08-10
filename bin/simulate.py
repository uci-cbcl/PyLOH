#!/usr/bin/env python

import sys
import numpy as np
import pysam
import time

centromere_starts = {}
centromere_ends = {}
centromere_starts['chr1'] = 121100000
centromere_ends['chr1'] = 128000000

chrom_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

chrom_lens = [247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
              158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
              114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
              63811651, 62435964, 46944323, 49691432]

reads_len = 100
eps = 1e-16
c_N = 2

phred = 'I'*reads_len

innormal_paternal = '../normal.paternal.fasta'
innormal_maternal = '../normal.maternal.fasta'
intumor_paternal = '../tumor.paternal.fasta'
intumor_maternal = '../tumor.maternal.fasta'


def main():
    inbed = sys.argv[1]
    outnormal = sys.argv[2]
    outtumor = sys.argv[3]
    tumor_frec = float(sys.argv[4])
    reads_num = int(sys.argv[5])
    
    inbed = open(inbed)
    outnormal = open(outnormal, 'w')
    outtumor = open(outtumor, 'w')
    
    normal_p_fasta = pysam.Fastafile(innormal_paternal)
    normal_m_fasta = pysam.Fastafile(innormal_maternal)
    tumor_p_fasta = pysam.Fastafile(intumor_paternal)
    tumor_m_fasta = pysam.Fastafile(intumor_maternal)
    
    chroms = []
    starts = []
    ends = []
    seg_lens = []
    p_copies = []
    m_copies = []
    c_Es = []
    for line in inbed:
        if line[0] == '#':
            continue
        
        fields = line.strip('\n').split('\t')
        chrom = fields[0]
        start, end, p_copy, m_copy = map(int, fields[1:5])
        
        chroms.append(chrom)
        starts.append(start)
        ends.append(end)
        seg_lens.append(end - start)
        p_copies.append(p_copy)
        m_copies.append(m_copy)
        
        c_T = p_copy + m_copy
        c_E = (1 - tumor_frec)*c_N + tumor_frec*c_T
        c_Es.append(c_E)
        
    seg_lens = np.array(seg_lens)
    c_Es = np.array(c_Es)
    
    seg_P = seg_lens*c_Es
    seg_P = seg_P/seg_P.sum()
    
    time_start = time.time()
    time_current = time.time()
    
    #simulate normal reads
    i = 0
    while i < reads_num:
        j = np.random.multinomial(1, seg_P, size=1)[0].tolist().index(1)
        chrom = chroms[j]
        chrom_idx = chrom_list.index(chrom)
        start = 0
        end = chrom_lens[chrom_idx]
        p_paternal = 0.5
        pos = np.random.random_integers(start, end - reads_len, 1)[0]
        
        if pos >= centromere_starts[chrom] and pos <= centromere_ends[chrom]:
            continue
        
        if np.random.rand() <= 0.5:
            seq = normal_p_fasta.fetch(chrom, pos, pos + reads_len)
            flag = 'paternal'
        else:
            seq = normal_m_fasta.fetch(chrom, pos, pos + reads_len)
            flag = 'maternal'
        
        outnormal.write('@{0}:{1}:{2}:{3}\n'.format(i + 1, flag, chrom, pos))
        outnormal.write('\n'.join([seq, '+', phred]) + '\n')
        
        i += 1
        
        if i % 100000 == 0:
            time_runing = time.time() - time_current
            time_current = time.time()
            print '{0} normal reads simulated...{1} sec'.format(i, time_runing)
            sys.stdout.flush()
    
    #simulate tumor reads
    for i in range(0, reads_num):
        j = np.random.multinomial(1, seg_P, size=1)[0].tolist().index(1)
        chrom = chroms[j]
        start = starts[j]
        end = ends[j]
        p_copy = p_copies[j]
        m_copy = m_copies[j]
        c_T = p_copy + m_copy
        p_tumor = tumor_frec*c_T/(tumor_frec*c_T + (1 - tumor_frec)*c_N)
        if p_copy == 0 and m_copy == 0:
            p_paternal = 0.5
        else:
            p_paternal = p_copy*1.0/(p_copy + m_copy)
                
        pos = np.random.random_integers(start, end - reads_len, 1)[0]
        
        if np.random.rand() <= p_tumor:
            if np.random.rand() <= p_paternal:
                seq = tumor_p_fasta.fetch(chrom, pos, pos + reads_len)
                flag = 'paternal'
            else:
                seq = tumor_m_fasta.fetch(chrom, pos, pos + reads_len)
                flag = 'maternal'
        else:
            if np.random.rand() <= 0.5:
                seq = normal_p_fasta.fetch(chrom, pos, pos + reads_len)
                flag = 'paternal'
            else:
                seq = normal_m_fasta.fetch(chrom, pos, pos + reads_len)
                flag = 'maternal'
        
        outtumor.write('@{0}:{1}:{2}:{3}\n'.format(i + 1, flag, chrom, pos))
        outtumor.write('\n'.join([seq, '+', phred]) + '\n')
        
        if i % 100000 == 0:
            time_runing = time.time() - time_current
            time_current = time.time()
            print '{0} tumor reads simulated...{1} sec'.format(i, time_runing)
            sys.stdout.flush()

    time_end = time.time()
    
    print 'Total time: {0}'.format(time_end - time_start)

    inbed.close()
    outnormal.close()
    outtumor.close()
    
if __name__ == '__main__':
    main()
