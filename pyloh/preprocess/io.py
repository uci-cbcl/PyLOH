'''
Created

@author: Andrew Roth

JointSNVMix-0.6.2
joint_snv_mix.pre_processing.bam_to_jcnt.BamToJcntConverter
joint_snv_mix.pre_processing.bam_to_jcnt.JointPileupIterator

================================================================================

Modified on 2013-07-20

@author: Yi Li
'''
import sys
from collections import Counter

import numpy as np
import pysam

ascii_offset = 33

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
            
            if ref_base == '':
                print 'Error: %s does not match the reference of the bam files' \
                % self.ref_genome_fasta.filename
                sys.exit(-1)
            
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


    
    
    
    
