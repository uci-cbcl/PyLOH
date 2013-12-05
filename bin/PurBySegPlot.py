#!/usr/bin/env python

import sys
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

def main():
    infile = open(sys.argv[1])
    PyLOH = float(sys.argv[2])
    ABSOLUTE = float(sys.argv[3])
    
    purity_lst = []
    copy_number_lst = []
    
    for line in infile:
        if line[0] == '#':
            continue
        
        fields = line.strip('\n').split('\t')
        LOH_status = fields[3]
        purity = float(fields[5])
        copy_number = int(fields[6])
        
        if LOH_status == 'TRUE':
            purity_lst.append(purity)
            copy_number_lst.append(copy_number)
            
    seg_num = len(purity_lst)
    
    plt.plot(range(seg_num), purity_lst, 'bs')
    plt.plot([0, seg_num-1], [PyLOH, PyLOH], '--', color = 'green', label = 'PyLOH')
    plt.plot([0, seg_num-1], [ABSOLUTE, ABSOLUTE], '-', color = 'red', label = 'ABSOLUTE')
    plt.xlim(0, seg_num-1)
    plt.ylim(0, 1)
    plt.xticks(sp.linspace(0, seg_num-1, seg_num), copy_number_lst)
    plt.yticks(sp.linspace(0, 1, 11), ['0%', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%', '100%'])
    plt.xlabel('Copy Number')
    plt.ylabel('Purity')
    plt.legend(loc=4,prop={'size':12})
    plt.show()
    
    infile.close()

if __name__ == '__main__':
    main()