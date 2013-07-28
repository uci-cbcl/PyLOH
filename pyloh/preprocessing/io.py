'''
Created on 2013-07-20

@author: Yi Li
'''
import numpy as np
from collections import Counter

from pyloh import constants

ascii_offset = 33

class Data:
    def __init__(self, segments=None, paired_counts=None):
        self.segments = Segments()
        self.paired_counts = []
        self.seg_num = 0
        
        if segments != None and paired_counts != None:
            self.segments = segments
            self.paired_counts = paired_counts
            self.seg_num = segments.num

    def read_data(self, filename_base):
        inseg_file_name = filename_base + '.segments'
        incounts_file_name = filename_base + '.counts'
            
        self.segments.read_segfile(inseg_file_name)
        self.seg_num = self.segments.num
        self._read_counts(incounts_file_name)
        
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

class Segments:
    def __init__(self):
        self.num = 0
        self.names = []
        self.chroms = []
        self.starts = []
        self.ends = []
        self.normal_reads_num = []
        self.tumor_reads_num = []
        
    def segmentation_by_chrom(self, normal_bam, tumor_bam):
        chrom_list = constants.CHROM_LIST
        chrom_lens = constants.CHROM_LENS
        chrom_start = constants.CHROM_START
        chrom_num = len(chrom_list)
        
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
        
    def segmentation_by_bed(self, normal_bam, tumor_bam, bed_file_name):
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

            normal_reads_num = normal_bam.count(bed_chroms[i], bed_starts[i], bed_ends[i])
            tumor_reads_num = tumor_bam.count(bed_chroms[i], bed_starts[i], bed_ends[i])
            
            self.names.append(seg_name)
            self.chroms.append(bed_chroms[i])
            self.starts.append(bed_starts[i])
            self.ends.append(bed_ends[i])
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            self.num = self.num + 1
    
    def read_segfile(self, inseg_file_name):
        infile = open(inseg_file_name)
        
        for line in infile:
            if line[0] == '#':
                continue
            
            fields = line.strip('\n').split('\t')
            seg_name, chrom = fields[0:2]
            start, end, normal_reads_num, tumor_reads_num = map(int, fields[2:6])
            
            self.names.append(seg_name)
            self.chroms.append(chrom)
            self.starts.append(start)
            self.ends.append(end)
            self.normal_reads_num.append(normal_reads_num)
            self.tumor_reads_num.append(tumor_reads_num)
            
            self.num = self.num + 1
        
        infile.close()
        
    def write_segfile(self, outseg_file_name):
        outfile = open(outseg_file_name, 'w')
        
        outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num', 'tumor_reads_num']) + '\n')
        
        for j in range(0, self.num):
            outfile.write('\t'.join(map(str, self[j])) + '\n')
        
        outfile.close()
    
    def __getitem__(self, i):
        "seg_name, chrom, start, end, normal_reads_num, tumor_reads_num"
        return (self.names[i], self.chroms[i], self.starts[i], self.ends[i], self.normal_reads_num[i], self.tumor_reads_num[i])
        
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
            ref_base = self.ref_genome_fasta.fetch(self.chrom, pos, pos + 1)
            
            paired_counts = self._get_paired_counts(normal_column, tumor_column, pos, ref_base)
            
            if paired_counts == None:
                normal_column, tumor_column = self.paired_pileup_iter.next()
                continue
            else:
                return paired_counts

    def _get_paired_counts(self, normal_column, tumor_column, pos, ref_base):
        normal_bases = self._parse_pileup_column(normal_column)
        tumor_bases = self._parse_pileup_column(tumor_column)
        
        tumor_non_ref_base, tumor_counts = self._get_counts(ref_base, tumor_bases)        
        normal_non_ref_base, normal_counts = self._get_counts(ref_base, normal_bases, non_ref_base=tumor_non_ref_base)        
    
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
            
            base = read.alignment.seq[qpos]
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
                raise Exception("Error in joint pileup iterator.")

#=======================================================================================================================
# Functions
#=======================================================================================================================

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
    
    
    
    