'''
Created on 2013-07-20

@author: Yi Li
'''
import sys
import time
from collections import Counter

import numpy as np
import pysam

from pyloh import constants
from pyloh.preprocess.data import Data, Segments
from pyloh.preprocess.utils import *

ascii_offset = 33

#JointSNVMix
def preprocess(args):        
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

#JointSNVMix
class PairedCountsIterator:
    def __init__(self, paired_pileup_iter, ref_genome_fasta, chrom, min_depth=20, min_bqual=10, min_mqual=10):
        self.paired_pileup_iter = paired_pileup_iter
        self.ref_genome_fasta = ref_genome_fasta
        self.chrom = chrom
        self.min_depth = min_depth
        self.min_bqual = min_bqual
        self.min_mqual = min_mqual
        
    def __iter__(self):
        return self
    
    def next(self):
        normal_column, tumor_column = self.paired_pileup_iter.next()
        
        while True:
            if normal_column.n < self.min_depth or tumor_column.n < self.min_depth:
                normal_column, tumor_column = self.paired_pileup_iter.next()
                continue
            
            pos = normal_column.pos
            ref_base = self.ref_genome_fasta.fetch(self.chrom, pos, pos + 1).upper()
            
            paired_counts = self._get_paired_counts(normal_column, tumor_column, pos, ref_base)
            
            if paired_counts == None:
                normal_column, tumor_column = self.paired_pileup_iter.next()
                continue
            else:
                return paired_counts

    def _get_paired_counts(self, normal_column, tumor_column, pos, ref_base):
        normal_bases = self._parse_pileup_column(normal_column)
        tumor_bases = self._parse_pileup_column(tumor_column)
        
        normal_non_ref_base, normal_counts = self._get_counts(ref_base, normal_bases) 
        tumor_non_ref_base, tumor_counts = self._get_counts(ref_base, tumor_bases)        

        # Check again for lines below read depth. The first check above speeds things up, though redundant.
        normal_depth = normal_counts[0] + normal_counts[1]
        tumor_depth = tumor_counts[0] + tumor_counts[1]
    
        if normal_depth < self.min_depth or tumor_depth < self.min_depth:
            return None
    
        # Shift index to one based position.
        one_based_pos = pos + 1
        
        paired_counts = []
        paired_counts.extend(normal_counts)
        paired_counts.extend(tumor_counts)
        
        return paired_counts
               
    def _parse_pileup_column(self, pileup_column):
        bases = []
        
        for read in pileup_column.pileups:
            if read.is_del:
                continue
            
            qpos = read.qpos
                 
            mqual = read.alignment.mapq
            
            if mqual < self.min_mqual:
                continue            
            
            bqual = ord(read.alignment.qual[qpos]) - ascii_offset            
            
            if bqual < self.min_bqual:
                continue
            
            base = read.alignment.seq[qpos].upper()
            bases.append(base)
        
        return bases

    def _get_counts(self, ref_base, bases, non_ref_base=None):
        counter = Counter(bases)
        
        non_ref_base, counts = self._parse_counts(ref_base, counter, non_ref_base)
        
        return non_ref_base, counts
    
    def _parse_counts(self, ref_base, counter, non_ref_base=None):
        ref_counts = counter[ref_base]
        
        del counter[ref_base]
        del counter['N']
        
        # Check if there is any non-ref bases.
        if non_ref_base is not None:
            non_ref_counts = counter[non_ref_base]
        else:
            if len(counter) > 0:
                non_ref_base, non_ref_counts = counter.most_common(1)[0]
            else:
                non_ref_base = 'N'
                non_ref_counts = 0
        
        counts = (ref_counts, non_ref_counts)
        
        return non_ref_base, counts

#JointSNVMix
class PairedPileupIterator:
    def __init__(self, normal_iter, tumor_iter, segment_start, segment_end):
        self.normal_iter = normal_iter
        self.tumor_iter = tumor_iter
        self.segment_start = segment_start
        self.segment_end = segment_end
    
    def __iter__(self):
        return self
    
    def next(self):
        normal_column = self.normal_iter.next()
        tumor_column = self.tumor_iter.next()
        
        while True:                    
            normal_pos = normal_column.pos
            tumor_pos = tumor_column.pos
                  
            if normal_pos == tumor_pos:
                if normal_pos >= self.segment_start and normal_pos <= self.segment_end:
                    return normal_column, tumor_column
                else:
                    normal_column = self.normal_iter.next()
                    tumor_column = self.tumor_iter.next()
            elif normal_pos < tumor_pos:
                normal_column = self.normal_iter.next()
            elif normal_pos > tumor_pos:
                tumor_column = self.tumor_iter.next()
            else:
                raise Exception("Error in paired pileup iterator.")


    
    
    
    
