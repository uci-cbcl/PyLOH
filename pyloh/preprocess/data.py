'''
Created on 2013-08-13

@author: Yi Li
'''
import sys
import pickle
import numpy as np

from pyloh import constants
from pyloh.preprocess.utils import *

class Data:
    def __init__(self, segments=None, paired_counts=None):
        self.segments = Segments()
        self.paired_counts = []
        self.sites_num = []
        self.seg_num = 0
        
        if segments != None and paired_counts != None:
            self.segments = segments
            self.paired_counts = paired_counts
            self.sites_num = self._get_sites_num(paired_counts)
            self.seg_num = segments.num
        
    def _get_sites_num(self, paired_counts):
        sites_num = []
        seg_num = len(paired_counts)
        
        for j in range(0, seg_num):
            sites_num.append(paired_counts[j].shape[0])
            
        sites_num = np.array(sites_num)
        
        return sites_num

    def read_data(self, filename_base):
        inseg_file_name = filename_base + '.PyLOH.segments'
        incounts_file_name = filename_base + '.PyLOH.counts'
            
        self.segments.read_segfile(inseg_file_name)
        self.seg_num = self.segments.num
        self._read_counts(incounts_file_name)
        self.sites_num = self._get_sites_num(self.paired_counts)
        
    def _read_counts(self, incounts_file_name):
        infile = open(incounts_file_name)
        
        for j in range(0, self.seg_num):
            self.paired_counts.append([])
            
        for line in infile:
            if line[0] == '#':
                continue
            
            idx, a_N, b_N, a_T, b_T = map(int, line.strip('\n').split('\t'))
            
            self.paired_counts[idx].append([a_N, b_N, a_T, b_T])
            
        for j in range(0, self.seg_num):
            self.paired_counts[j] = np.array(self.paired_counts[j])
        
        infile.close()

    def write_data(self, filename_base):
        outseg_file_name = filename_base + '.PyLOH.segments'
        outcounts_file_name = filename_base + '.PyLOH.counts'
        
        self.segments.write_segfile(outseg_file_name)
        self._write_counts(outcounts_file_name)
        
    def _write_counts(self, outcounts_file_name):
        outfile = open(outcounts_file_name, 'w')
        
        outfile.write('\t'.join(['#seg_idx', 'normal_A', 'normal_B', 'tumor_A', 'tumor_B']) + '\n')
        
        for j in range(0, self.seg_num):
            for i in range(0, self.paired_counts[j].shape[0]):
                outfile.write(str(j) + '\t' + '\t'.join(map(str, self.paired_counts[j][i])) + '\n')
        
        outfile.close()
        
    def tumor_LOH_test(self, WES_flag):
        self.segments.tumor_LOH_test(self.paired_counts, WES_flag)

class Segments:
    def __init__(self):
        self.num = 0
        self.Lambda_S = 0
        self.names = []
        self.chroms = []
        self.starts = []
        self.ends = []
        self.normal_reads_num = []
        self.tumor_reads_num = []
        self.LOH_frec = []
        self.LOH_status = []
        self.log2_ratio = []
        
    def segmentation_by_chrom(self, normal_bam, tumor_bam):
        chrom_ID_list = constants.CHROM_ID_LIST
        chrom_start = constants.CHROM_START
        chrom_num = len(chrom_ID_list)
        
        sam_SQ = normal_bam.header['SQ']
        sam_chrom_format = get_chrom_format(map(lambda x:x['SN'], sam_SQ))
        chrom_lens = self._get_chrom_lens(chrom_ID_list, sam_SQ)
        
        for i in range(0, chrom_num):
            chrom_name = chrom_ID_to_name(chrom_ID_list[i], sam_chrom_format)
            seg_name = self._get_segment_name(chrom_name, chrom_start, chrom_lens[i])
            normal_reads_num = normal_bam.count(chrom_name, chrom_start, chrom_lens[i])
            tumor_reads_num = tumor_bam.count(chrom_name, chrom_start, chrom_lens[i])
            
            self.names.append(seg_name)
            self.chroms.append(chrom_name)
            self.starts.append(chrom_start)
            self.ends.append(chrom_lens[i])
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.log2_ratio.append(np.log2(1.0*tumor_reads_num/normal_reads_num))
        
        self.num = chrom_num
        self._init_LOH_status()
        
    def segmentation_by_bed(self, normal_bam, tumor_bam, bed_file_name):
        chrom_ID_list = constants.CHROM_ID_LIST
        chrom_start = constants.CHROM_START
        
        sam_SQ = normal_bam.header['SQ']
        sam_chrom_format = get_chrom_format(map(lambda x:x['SN'], sam_SQ))
        chrom_lens = self._get_chrom_lens(chrom_ID_list, sam_SQ)
        
        bed_chroms, bed_starts, bed_ends = BEDParser(bed_file_name)
        bed_chrom_format = get_chrom_format(bed_chroms)
        bed_num = len(bed_chroms)
        
        for i in range(0, bed_num):
            chrom_ID = chrom_name_to_ID(bed_chroms[i])
            chrom_name = chrom_ID_to_name(chrom_ID, sam_chrom_format)
            seg_name = self._get_segment_name(chrom_name, bed_starts[i], bed_ends[i])
            
            if chrom_ID not in chrom_ID_list:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue
            
            chrom_idx = chrom_ID_list.index(chrom_ID)
            
            if bed_starts[i] < chrom_start or bed_ends[i] > chrom_lens[chrom_idx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue

            normal_reads_num = normal_bam.count(chrom_name, bed_starts[i], bed_ends[i])
            tumor_reads_num = tumor_bam.count(chrom_name, bed_starts[i], bed_ends[i])
            
            self.names.append(seg_name)
            self.chroms.append(chrom_name)
            self.starts.append(bed_starts[i])
            self.ends.append(bed_ends[i])
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.log2_ratio.append(np.log2(1.0*tumor_reads_num/normal_reads_num))
            self.num = self.num + 1
            
        self._init_LOH_status()
    
    def _get_chrom_lens(self, chrom_ID_list, sam_SQ):
        chrom_lens = []
        
        for i in range(0, len(chrom_ID_list)):
            chrom_ID = chrom_ID_list[i]
            
            for j in range(0, len(sam_SQ)):
                if chrom_ID == chrom_name_to_ID(sam_SQ[j]['SN']):
                    chrom_lens.append(int(sam_SQ[j]['LN']))
                    break
        
        return chrom_lens
        
    def tumor_LOH_test(self, paired_counts, WES_flag):
        for j in range(0, self.num):
            LOH_frec, LOH_status = tumor_LOH_test(paired_counts[j], WES_flag)
            self.LOH_frec[j] = LOH_frec
            self.LOH_status[j] = LOH_status
            
    def compute_Lambda_S(self):
        reads_depth_ratio = []
        
        for j in range(0, self.num):
            if self.LOH_status[j] == 'FALSE':
                ratio = 1.0*self.tumor_reads_num[j]/self.normal_reads_num[j]
                reads_depth_ratio.append(ratio)
                
        reads_depth_ratio = np.array(reads_depth_ratio)
        reads_depth_ratio = remove_outliers(reads_depth_ratio)
        
        if reads_depth_ratio.shape[0] != 0:
            self.Lambda_S = reads_depth_ratio.mean()
        else:
            print 'Error: No diploid segments found, existing...'
            sys.exit(1)
        
    def read_segfile(self, inseg_file_name):
        infile = open(inseg_file_name)
        
        for line in infile:
            if line[0] == '#':
                continue
            
            fields = line.strip('\n').split('\t')
            seg_name, chrom = fields[0:2]
            start, end, normal_reads_num, tumor_reads_num = map(int, fields[2:6])
            LOH_frec = float(fields[6])
            LOH_status = fields[7]
            log2_ratio = float(fields[8])
            
            self.names.append(seg_name)
            self.chroms.append(chrom)
            self.starts.append(start)
            self.ends.append(end)
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.LOH_frec.append(LOH_frec)
            self.LOH_status.append(LOH_status)
            self.log2_ratio.append(log2_ratio)
            
            self.num = self.num + 1
        
        infile.close()
        
    def write_segfile(self, outseg_file_name):
        outfile = open(outseg_file_name, 'w')
        
        outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num',
                                 'tumor_reads_num', 'LOH_frec', 'LOH_status', 'log2_ratio']) + '\n')
        
        for j in range(0, self.num):
            outfile.write('\t'.join(map(str, self[j])) + '\n')
        
        outfile.close()
    
    def _init_LOH_status(self):
        self.LOH_frec = ['NONE' for j in range(0, self.num)]
        self.LOH_status = ['NONE' for j in range(0, self.num)]
    
    def __getitem__(self, i):
        "seg_name, chrom, start, end, normal_reads_num, tumor_reads_num, LOH_frec, LOH_status, log2_ratio"
        return (self.names[i], self.chroms[i], self.starts[i], self.ends[i],
                self.normal_reads_num[i], self.tumor_reads_num[i],
                self.LOH_frec[i], self.LOH_status[i], self.log2_ratio[i])
        
    def _get_segment_name(self, chrom, start, end):
        return '_'.join([chrom, 'start', str(start), 'end', str(end)])
    
class BAFHeatMap:
    def __init__(self, BAF_counts=None):
        self.BAF_counts = []
        self.seg_num = 0
        
        if BAF_counts != None:
            self.BAF_counts = BAF_counts
            self.seg_num = len(BAF_counts)
    
    def write_heatmap(self, filename_base):
        outheatmap_file_name = filename_base + '.PyLOH.heatmap.pkl'
        outfile = open(outheatmap_file_name, 'wb')
        
        pickle.dump(self.BAF_counts, outfile)
        
        outfile.close()
        
    def read_heatmap(self, filename_base):
        inheatmap_file_name = filename_base + '.PyLOH.heatmap.pkl'
        infile = open(inheatmap_file_name, 'rb')
        
        self.BAF_counts = pickle.load(infile)
        self.seg_num = len(self.BAF_counts)
        
        infile.close()
        
    def get_color_max(self):
        BAF_counts_min = constants.BAF_COUNTS_MIN
        BAF_counts_max = constants.BAF_COUNTS_MAX
        
        color_max = 0
        for j in range(0, self.seg_num):
            BAF_counts_j_sub = self.BAF_counts[j][BAF_counts_min:BAF_counts_max, BAF_counts_min:BAF_counts_max]
            color_max_j = BAF_counts_j_sub.max()
            color_max = max(color_max, color_max_j)
            
        self.color_max = color_max
        
        
        
        