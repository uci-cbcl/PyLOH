#!/usr/bin/env python

#=======================================================================================================================
# PyLOH
# Author : Yi Li
# E-mail : yil8@uci.edu
#=======================================================================================================================
import argparse

from pyloh.preprocessing.io import preprocess

parser = argparse.ArgumentParser(prog='PyLOH')
subparsers = parser.add_subparsers()

#===============================================================================
# Add bam2data sub-command
#===============================================================================
parser_jcnt = subparsers.add_parser('preprocess',
                                    help='Preprocess paired normal and tumor BAM files')

parser_jcnt.add_argument('reference_genome_file_name',
                          help='''Reference genome fasta file.''')

parser_jcnt.add_argument('normal_bam_file_name',
                          help='''Normal BAM file.''')

parser_jcnt.add_argument('tumor_bam_file_name',
                          help='''Tumor BAM file.''')

parser_jcnt.add_argument('data_file_basename',
                          help='Base name of preprocessed files to be created.')

parser_jcnt.add_argument('--segments_bed_file_name', default=None, type=str,
                          help='''Bed file for segments. If not provided, use autosomes as segments. Default is None.''')

parser_jcnt.add_argument('--min_depth', default=20, type=int,
                          help='''Minimum depth of coverage in both tumor and normal sample required to use a site in
                          the analysis. Default is 20.''')

parser_jcnt.add_argument('--min_base_qual', default=10, type=int,
                          help='''Remove bases with base quality lower than this. Default is 10.''')

parser_jcnt.add_argument('--min_map_qual', default=10, type=int,
                          help='''Remove bases with mapping quality lower than this. Default is 10.''')

parser_jcnt.set_defaults(func=preprocess)





#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)