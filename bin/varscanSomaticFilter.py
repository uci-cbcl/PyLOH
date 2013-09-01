#!/usr/bin/env python

import sys

def main():
    infile = open(sys.argv[1])
    outfile = open(sys.argv[2], 'w')
    depth_thred = int(sys.argv[3])
    
    for line in infile:
        if line[0:5] == 'chrom':
            continue
        
        fields = line.split('\t')
        a_T = int(fields[8])
        b_T = int(fields[9])
        d_T = a_T + b_T
        somatic_status = fields[12]
        
        if d_T >= depth_thred and somatic_status == 'Somatic':
            outfile.write(line)
    
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main()