'''
Created on 2013-08-14

@author: Yi Li
'''
import sys
import time

import numpy as np
import pysam

from pyloh import constants
from pyloh.preprocess.data import Data, Segments
from pyloh.preprocess.io import PairedCountsIterator, PairedPileupIterator
from pyloh.preprocess.utils import *

#JointSNVMix
def run_preprocess(args):        
    normal_bam = pysam.Samfile(args.normal_bam_file_name, 'rb')
    tumor_bam = pysam.Samfile(args.tumor_bam_file_name, 'rb')
    
    ref_genome_fasta = pysam.Fastafile(args.reference_genome_file_name)
    
    segments = Segments()
    
    if args.segments_bed_file_name == None:
        print 'Loading segments by 22 autosomes...'
        sys.stdout.flush()
        segments.segmentation_by_chrom(normal_bam, tumor_bam)
    else:
        print 'Loading segments by {0}...'.format(args.segments_bed_file_name)
        sys.stdout.flush()
        segments.segmentation_by_bed(normal_bam, tumor_bam, args.segments_bed_file_name)
    
    time_start = time.time()           
    
    converter = BamToDataConverter(
                                   normal_bam,
                                   tumor_bam,
                                   ref_genome_fasta,
                                   args.data_file_basename,
                                   segments,
                                   min_depth=args.min_depth,
                                   min_bqual=args.min_base_qual,
                                   min_mqual=args.min_map_qual,
                                   )
    
    converter.convert()
    
    time_end = time.time()
    
    print 'Run time: {0:.2f} seconds'.format(time_end - time_start)
    sys.stdout.flush()

#JointSNVMix
class BamToDataConverter:
    def __init__(self, normal_bam, tumor_bam, ref_genome_fasta, data_file_basename,
                 segments, min_depth=20, min_bqual=10, min_mqual=10):
        self.normal_bam = normal_bam
        self.tumor_bam = tumor_bam
        self.ref_genome_fasta = ref_genome_fasta
        self.data_file_basename = data_file_basename
        
        self.segments = segments
        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        
        self.buffer = 100000
        
    def convert(self):
        seg_num = self.segments.num
        paired_counts = []
        
        for j in range(0, seg_num):
            print 'Preprocessing segment {0}...'.format(self.segments[j][0])
            sys.stdout.flush()
            
            chrom = self.segments[j][1]
            start = self.segments[j][2]
            end = self.segments[j][3]
            
            normal_pileup_iter = self.normal_bam.pileup(chrom, start, end)
            tumor_pileup_iter = self.tumor_bam.pileup(chrom, start, end)
            
            paired_pileup_iter = PairedPileupIterator(normal_pileup_iter, tumor_pileup_iter, start, end)
            paired_counts_iter = PairedCountsIterator(paired_pileup_iter, self.ref_genome_fasta, chrom,
                                                      self.min_depth, self.min_bqual, self.min_mqual)
            paired_counts_j = self._convert_by_segments(paired_counts_iter)
            
            paired_counts.append(paired_counts_j)
        
        data = Data(self.segments, paired_counts)
        data.tumor_LOH_test()
        data.write_data(self.data_file_basename)
        
    def _convert_by_segments(self, paired_counts_iter):
        paired_counts_j = np.array([[], [], [], []], dtype=int).transpose()
        buffer_counts = []
        i = 0
        
        for counts in paired_counts_iter:
            buffer_counts.append(counts)
            i = i + 1
            
            if i < self.buffer:
                continue
            
            buffer_counts = np.array(buffer_counts)
            buffer_counts_filtered = normal_heterozygous_filter(buffer_counts)
            if buffer_counts_filtered.shape[0] != 0:
                paired_counts_j = np.vstack((paired_counts_j, buffer_counts_filtered))
            
            buffer_counts = []
            i = 0
        
        buffer_counts = np.array(buffer_counts)
        buffer_counts_filtered = normal_heterozygous_filter(buffer_counts)
        if buffer_counts_filtered.shape[0] != 0:
            paired_counts_j = np.vstack((paired_counts_j, buffer_counts_filtered))
        
        return paired_counts_j