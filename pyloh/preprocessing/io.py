'''
Created on 2013-07-20

@author: Yi Li
'''
import numpy as np

from pyloh import constants



class Segments:
    def __init__(self):
        self.num = 0
        self.names = []
        self.chroms = []
        self.starts = []
        self.ends = []
        
    def segmentation_by_chrom(self):
        chrom_list = constants.CHROM_LIST
        chrom_lens = constants.CHROM_LENS
        chrom_start = constants.CHROM_START
        chrom_num = len(chrom_list)
        
        for i in range(0, chrom_num):
            seg_name = self._get_segment_name(chrom_list[i], chrom_start, chrom_lens[i])
            
            self.names.append(seg_name)
            self.chroms.append(chrom_list[i])
            self.starts.append(chrom_start)
            self.ends.append(chrom_lens[i])
        
        self.num = chrom_num
        
    def segmentation_by_bed(self, bed_file_name):
        chrom_list = constants.CHROM_LIST
        chrom_lens = constants.CHROM_LENS
        chrom_start = constants.CHROM_START
        
        bed_chroms, bed_starts, bed_ends = BEDParser(bed_file_name)
        bed_num = len(bed_chroms)
        
        for i in range(0, bed_num):
            seg_name = self._get_segment_name(bed_chroms[i], bed_starts[i], bed_ends[i])
            
            if bed_chroms[i] not in chrom_list:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                continue
            
            if bed_starts[i] < chrom_start or bed_ends[i] > chrom_lens[i]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                continue
            
            self.names.append(seg_name)
            self.chroms.append(bed_chroms[i])
            self.starts.append(bed_starts[i])
            self.ends.append(bed_ends[i])
            self.num = self.num + 1
    
    def __getitem__(self, i):
        return (self.names[i], self.chroms[i], self.starts[i], self.ends[i])
        
    def _get_segment_name(self, chrom, start, end):
        return '_'.join([chrom, 'start', str(start), 'end', str(end)])
        


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
    
    
    
    