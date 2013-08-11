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
from pyloh.preprocessing.utils import *

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
        inseg_file_name = filename_base + '.segments'
        incounts_file_name = filename_base + '.counts'
            
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
        outseg_file_name = filename_base + '.segments'
        outcounts_file_name = filename_base + '.counts'
        
        self.segments.write_segfile(outseg_file_name)
        self._write_counts(outcounts_file_name)
        
    def _write_counts(self, outcounts_file_name):
        outfile = open(outcounts_file_name, 'w')
        
        outfile.write('\t'.join(['#seg_idx', 'normal_A', 'normal_B', 'tumor_A', 'tumor_B']) + '\n')
        
        for j in range(0, self.seg_num):
            for i in range(0, self.paired_counts[j].shape[0]):
                outfile.write(str(j) + '\t' + '\t'.join(map(str, self.paired_counts[j][i])) + '\n')
        
        outfile.close()
        
    def tumor_LOH_test(self):
        self.segments.tumor_LOH_test(self.paired_counts)

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
        
    def segmentation_by_chrom(self, normal_bam, tumor_bam):
        chrom_list = constants.CHROM_LIST
        chrom_start = constants.CHROM_START
        chrom_num = len(chrom_list)
        
        header_SQ = normal_bam.header['SQ']
        chrom_lens = self._get_chrom_lens(chrom_list, header_SQ)
        
        for i in range(0, chrom_num):
            seg_name = self._get_segment_name(chrom_list[i], chrom_start, chrom_lens[i])
            normal_reads_num = normal_bam.count(chrom_list[i], chrom_start, chrom_lens[i])
            tumor_reads_num = tumor_bam.count(chrom_list[i], chrom_start, chrom_lens[i])
            
            self.names.append(seg_name)
            self.chroms.append(chrom_list[i])
            self.starts.append(chrom_start)
            self.ends.append(chrom_lens[i])
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
        
        self.num = chrom_num
        self._init_LOH_status()
        
    def segmentation_by_bed(self, normal_bam, tumor_bam, bed_file_name):
        chrom_list = constants.CHROM_LIST
        chrom_start = constants.CHROM_START
        
        header_SQ = normal_bam.header['SQ']
        chrom_lens = self._get_chrom_lens(chrom_list, header_SQ)
        
        bed_chroms, bed_starts, bed_ends = BEDParser(bed_file_name)
        bed_num = len(bed_chroms)
        
        for i in range(0, bed_num):
            seg_name = self._get_segment_name(bed_chroms[i], bed_starts[i], bed_ends[i])
            
            if bed_chroms[i] not in chrom_list:
                print 'Chromsome {0} not found, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue
            
            chrom_idx = chrom_list.index(bed_chroms[i])
            
            if bed_starts[i] < chrom_start or bed_ends[i] > chrom_lens[chrom_idx]:
                print 'Out of range chromsome {0}, segment {1} excluded...'.format(bed_chroms[i], seg_name)
                sys.stdout.flush()
                continue

            normal_reads_num = normal_bam.count(bed_chroms[i], bed_starts[i], bed_ends[i])
            tumor_reads_num = tumor_bam.count(bed_chroms[i], bed_starts[i], bed_ends[i])
            
            self.names.append(seg_name)
            self.chroms.append(bed_chroms[i])
            self.starts.append(bed_starts[i])
            self.ends.append(bed_ends[i])
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.num = self.num + 1
            
        self._init_LOH_status()
    
    def _get_chrom_lens(self, chrom_list, header_SQ):
        chrom_lens = []
        
        for i in range(0, len(chrom_list)):
            chrom = chrom_list[i]
            
            for j in range(0, len(header_SQ)):
                if chrom == header_SQ[j]['SN']:
                    chrom_lens.append(int(header_SQ[j]['LN']))
                    break
        
        return chrom_lens
        
    def tumor_LOH_test(self, paired_counts):
        for j in range(0, self.num):
            LOH_frec, LOH_status = tumor_LOH_test(paired_counts[j])
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
            
            self.names.append(seg_name)
            self.chroms.append(chrom)
            self.starts.append(start)
            self.ends.append(end)
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.LOH_frec.append(LOH_frec)
            self.LOH_status.append(LOH_status)
            
            self.num = self.num + 1
        
        infile.close()
        
    def write_segfile(self, outseg_file_name):
        outfile = open(outseg_file_name, 'w')
        
        outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num',
                                 'tumor_reads_num', 'LOH_frec', 'LOH_status']) + '\n')
        
        for j in range(0, self.num):
            outfile.write('\t'.join(map(str, self[j])) + '\n')
        
        outfile.close()
    
    def _init_LOH_status(self):
        self.LOH_frec = ['NONE' for j in range(0, self.num)]
        self.LOH_status = ['NONE' for j in range(0, self.num)]
    
    def __getitem__(self, i):
        "seg_name, chrom, start, end, normal_reads_num, tumor_reads_num, LOH_frec, LOH_status"
        return (self.names[i], self.chroms[i], self.starts[i], self.ends[i],
                self.normal_reads_num[i], self.tumor_reads_num[i],
                self.LOH_frec[i], self.LOH_status[i])
        
    def _get_segment_name(self, chrom, start, end):
        return '_'.join([chrom, 'start', str(start), 'end', str(end)])

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


    
    
    
    
