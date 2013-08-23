#!/usr/bin/env python

#=======================================================================================================================
# PyLOH
# Author : Yi Li
# E-mail : yil8@uci.edu
# Notes  : Some design patterns are based on JointSNVMix (http://code.google.com/p/joint-snv-mix/)
#=======================================================================================================================
import argparse

from pyloh.preprocess.run_preprocess import run_preprocess
from pyloh.model.run_model import run_poisson_model
from pyloh.postprocess.plot import plot_BAF_heatmap


parser = argparse.ArgumentParser(prog='PyLOH')
subparsers = parser.add_subparsers()

#===============================================================================
# Add preprocess sub-command
#===============================================================================
parser_preprocess = subparsers.add_parser('preprocess',
                                    help='''Preprocess paired normal and tumor BAM files''')

parser_preprocess.add_argument('reference_genome_file_name',
                          help='''Reference genome fasta file.''')

parser_preprocess.add_argument('normal_bam_file_name',
                          help='''Normal BAM file.''')

parser_preprocess.add_argument('tumor_bam_file_name',
                          help='''Tumor BAM file.''')

parser_preprocess.add_argument('data_file_basename',
                          help='''Base name of preprocessed files to be created.''')

parser_preprocess.add_argument('--segments_bed_file_name', default=None, type=str,
                          help='''Bed file for segments. If not provided,
                            use autosomes as segments. Default is None.''')

parser_preprocess.add_argument('--min_depth', default=20, type=int,
                          help='''Minimum depth of coverage in both tumor and normal
                          sample required to use a site in the analysis. Default is 20.''')

parser_preprocess.add_argument('--min_base_qual', default=10, type=int,
                          help='''Remove bases with base quality lower than this. Default is 10.''')

parser_preprocess.add_argument('--min_map_qual', default=10, type=int,
                          help='''Remove bases with mapping quality lower than this. Default is 10.''')

parser_preprocess.add_argument('--process_num', default=1, type=int,
                          help='''Number of processes to launch for preprocessing. Default is 1.''')

parser_preprocess.set_defaults(func=run_preprocess)

#===============================================================================
# Add run_model sub-command
#===============================================================================
parser_run_model = subparsers.add_parser('run_model',
                                      help='''Run a Poisson model based analysis. Requires preprocessed counts
                                      file and segments file that have been created.''')

parser_run_model.add_argument('data_file_basename',
                            help='Base name of preprocessed files created.')

parser_run_model.add_argument('--priors_file_name', default=None, type=str,
                             help='''File containing prior distribution to use for training. If not provided,
                                use uniform prior. Default is None.''')

parser_run_model.add_argument('--max_iters', default=100, type=int,
                          help='''Maximum number of iterations to used for training model. Default 100''')

parser_run_model.add_argument('--stop_value', default=1e-10, type=float,
                          help='''Stop value for EM training. Once the change in log-likelihood function is below
                          this value training will end. Defaul 1e-10''')

parser_run_model.set_defaults(func=run_poisson_model)

#===============================================================================
# Add BAF_heatmap sub-command
#===============================================================================
parser_BAF_heatmap = subparsers.add_parser('BAF_heatmap',
                                      help='''Plot the BAF heat map for each segment. Requires preprocessed
                                      heatmap file that have been created.''')

parser_BAF_heatmap.add_argument('data_file_basename',
                            help='''Base name of preprocessed files created.''')

parser_BAF_heatmap.set_defaults(func=plot_BAF_heatmap)

#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)