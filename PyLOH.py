#!/usr/bin/env python

#=======================================================================================================================
# Created on 2013-07-28
# @author: Yi Li
#
# PyLOH
# Copyright (c) 2013 Yi Li <yil8@uci.edu>
#
# This code is free software; you can redistribute it and/or modify it
# under the terms of GNU GPL v2.0 (see the file LICENSE included with the distribution).

# Some features are built on top of JointSNVMix-0.6.2 (http://code.google.com/p/joint-snv-mix/).
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

parser_preprocess.add_argument('reference_genome',
                          help='''FASTA file for reference genome.''')

parser_preprocess.add_argument('normal_bam',
                          help='''BAM file for normal sample.''')

parser_preprocess.add_argument('tumor_bam',
                          help='''BAM file for tumor sample.''')

parser_preprocess.add_argument('filename_base',
                          help='''Base name of preprocessed files to be created.''')

parser_preprocess.add_argument('--segments_bed', default=None, type=str,
                          help='''BED file for segments. If not provided,
                            use autosomes as segments. Default is None.''')

parser_preprocess.add_argument('--min_depth', default=20, type=int,
                          help='''Minimum reads depth required for both normal and tumor samples. 
                          Default is 20.''')

parser_preprocess.add_argument('--min_base_qual', default=10, type=int,
                          help='''Minimum base quality required. Default is 10.''')

parser_preprocess.add_argument('--min_map_qual', default=10, type=int,
                          help='''Minimum mapping quality required. Default is 10.''')

parser_preprocess.add_argument('--process_num', default=1, type=int,
                          help='''Number of processes to launch for preprocessing. Default is 1.''')

parser_preprocess.set_defaults(func=run_preprocess)

#===============================================================================
# Add run_model sub-command
#===============================================================================
parser_run_model = subparsers.add_parser('run_model',
                                      help='''Run a probabilistic model based analysis. Requires preprocessed counts
                                      file and segments file that have been created.''')

parser_run_model.add_argument('filename_base',
                            help='Base name of preprocessed files created.')

parser_run_model.add_argument('--allele_number_max', default=2, type=int,
                            help='''Maximum copy number of each allele allows to take. Default is 2.''')

parser_run_model.add_argument('--priors_file_name', default=None, type=str,
                             help='''File of the prior distribution. If not provided,
                                use uniform prior. Default is None.''')

parser_run_model.add_argument('--max_iters', default=100, type=int,
                          help='''Maximum number of iterations for training. Default is 100.''')

parser_run_model.add_argument('--stop_value', default=1e-7, type=float,
                          help='''Stop value of the EM algorithm for training. If the change of log-likelihood is lower
                          than this value, stop training. Default is 1e-7.''')

parser_run_model.set_defaults(func=run_poisson_model)

#===============================================================================
# Add BAF_heatmap sub-command
#===============================================================================
parser_BAF_heatmap = subparsers.add_parser('BAF_heatmap',
                                      help='''Plot the BAF heat map for each segment. Requires preprocessed
                                      heatmap file that have been created.''')

parser_BAF_heatmap.add_argument('filename_base',
                            help='''Base name of preprocessed files created.''')

parser_BAF_heatmap.set_defaults(func=plot_BAF_heatmap)

#===============================================================================
# Run
#===============================================================================
args = parser.parse_args()

args.func(args)
