'''
Created on 2013-08-14

@author: Yi Li
'''
import sys
import time
from multiprocessing import Pool

import numpy as np
import pysam

from pyloh import constants
from pyloh.preprocess.data import Data, Segments, BAFHeatMap
from pyloh.preprocess.io import PairedCountsIterator, PairedPileupIterator
from pyloh.preprocess.utils import *

def run_preprocess(args):        
    normal_bam = pysam.Samfile(args.normal_bam_file_name, 'rb')
    tumor_bam = pysam.Samfile(args.tumor_bam_file_name, 'rb')
    
    segments = Segments()
    
    if args.segments_bed_file_name == None:
        print 'Loading segments by 22 autosomes...'
        sys.stdout.flush()
        segments.segmentation_by_chrom(normal_bam, tumor_bam)
    else:
        print 'Loading segments by {0}...'.format(args.segments_bed_file_name)
        sys.stdout.flush()
        segments.segmentation_by_bed(normal_bam, tumor_bam, args.segments_bed_file_name)
    
    normal_bam.close()
    tumor_bam.close()
    
    time_start = time.time()           
    
    converter = BamToDataConverter(
                                   args.normal_bam_file_name,
                                   args.tumor_bam_file_name,
                                   args.reference_genome_file_name,
                                   args.filename_base,
                                   segments,
                                   min_depth=args.min_depth,
                                   min_bqual=args.min_base_qual,
                                   min_mqual=args.min_map_qual,
                                   process_num = args.process_num
                                   )
    
    converter.convert()
    
    time_end = time.time()
    
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()

class BamToDataConverter:
    def __init__(self, normal_bam_file_name, tumor_bam_file_name,
                 reference_genome_file_name, filename_base,
                 segments, min_depth=20, min_bqual=10, min_mqual=10, process_num=1):
        self.normal_bam_file_name = normal_bam_file_name
        self.tumor_bam_file_name = tumor_bam_file_name
        self.reference_genome_file_name = reference_genome_file_name
        self.filename_base = filename_base
        
        self.segments = segments
        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        self.process_num = process_num
        
    def convert(self):
        seg_num = self.segments.num
        
        process_num = self.process_num
        
        if process_num > seg_num:
            process_num = seg_num
        
        pool = Pool(processes = process_num)
        
        args_list = []
        
        for j in range(0, seg_num):
            seg_name = self.segments[j][0]
            chrom = self.segments[j][1]
            start = self.segments[j][2]
            end = self.segments[j][3]
            
            args_tuple = (seg_name, chrom, start, end, self.normal_bam_file_name,
                          self.tumor_bam_file_name, self.reference_genome_file_name,
                          self.min_depth, self.min_bqual, self.min_mqual)
            
            args_list.append(args_tuple)
            
        counts_tuple_list = pool.map(process_by_segment, args_list)
        
        paired_counts = []
        BAF_counts = []
        
        for counts_tuple_j in counts_tuple_list:
            paired_counts_j, BAF_counts_j = counts_tuple_j
            paired_counts.append(paired_counts_j)
            BAF_counts.append(BAF_counts_j)
        
        BAF_heatmap = BAFHeatMap(BAF_counts)
        BAF_heatmap.write_heatmap(self.filename_base)
        
        data = Data(self.segments, paired_counts)
        data.tumor_LOH_test()
        data.write_data(self.filename_base)
    
#===============================================================================
# Function
#===============================================================================
def process_by_segment(args_tuple):
    seg_name, chrom, start, end, normal_bam_file_name, tumor_bam_file_name, \
    reference_genome_file_name, min_depth, min_bqual, min_mqual= args_tuple
    
    print 'Preprocessing segment {0}...'.format(seg_name)
    sys.stdout.flush()

    normal_bam = pysam.Samfile(normal_bam_file_name, 'rb')
    tumor_bam = pysam.Samfile(tumor_bam_file_name, 'rb')
    ref_genome_fasta = pysam.Fastafile(reference_genome_file_name)
    
    normal_pileup_iter = normal_bam.pileup(chrom, start, end)
    tumor_pileup_iter = tumor_bam.pileup(chrom, start, end)
    
    paired_pileup_iter = PairedPileupIterator(normal_pileup_iter, tumor_pileup_iter, start, end)
    paired_counts_iter = PairedCountsIterator(paired_pileup_iter, ref_genome_fasta, chrom,
                                              min_depth, min_bqual, min_mqual)
    
    paired_counts_j, BAF_counts_j = iterator_to_counts(paired_counts_iter)
    counts_tuple_j = (paired_counts_j, BAF_counts_j)
    
    normal_bam.close()
    tumor_bam.close()
    ref_genome_fasta.close()
    
    return counts_tuple_j

def iterator_to_counts(paired_counts_iter):
    buffer = 100000
    
    paired_counts_j = np.array([[], [], [], []], dtype=int).transpose()
    BAF_counts_j = np.zeros((100, 100))
    buffer_counts = []
    i = 0
        
    for counts in paired_counts_iter:
        buffer_counts.append(counts)
        i = i + 1
            
        if i < buffer:
            continue
            
        buffer_counts = np.array(buffer_counts)
        
        if buffer_counts.shape[0] != 0 :
            BAF_counts_buffer = get_BAF_counts(buffer_counts)
            BAF_counts_j += BAF_counts_buffer
        
        buffer_counts_filtered = normal_heterozygous_filter(buffer_counts)
        
        if buffer_counts_filtered.shape[0] != 0:
            paired_counts_j = np.vstack((paired_counts_j, buffer_counts_filtered))
            
        buffer_counts = []
        i = 0
        
    buffer_counts = np.array(buffer_counts)
    
    if buffer_counts.shape[0] != 0 :
        BAF_counts_buffer = get_BAF_counts(buffer_counts)
        BAF_counts_j += BAF_counts_buffer
    
    buffer_counts_filtered = normal_heterozygous_filter(buffer_counts)
    
    if buffer_counts_filtered.shape[0] != 0:
        paired_counts_j = np.vstack((paired_counts_j, buffer_counts_filtered))
        
    return (paired_counts_j, BAF_counts_j)
