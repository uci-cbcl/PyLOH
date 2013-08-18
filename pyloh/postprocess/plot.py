'''
Created on 2012-08-18

@author: Yi Li
'''
import os

import numpy as np
import scipy
from matplotlib import pyplot

from pyloh import constants
from pyloh.preprocess.data import Segments, BAFHeatMap

def plot_BAF_heatmap(args):
    BAF_heatmap = BAFHeatMap()
    BAF_heatmap.read_heatmap(args.data_file_basename)
    BAF_heatmap.get_color_max()
    
    segments = Segments()
    inseg_file_name = args.data_file_basename + '.PyLOH.segments'
    segments.read_segfile(inseg_file_name)
    seg_num = segments.num
    
    outheatmap_dir_name = args.data_file_basename + '.PyLOH.heatmap.plot'
    if os.path.exists(outheatmap_dir_name) == False:
        os.mkdir(outheatmap_dir_name)
    
    for j in range(0, seg_num):
        plot_BAF_heatmap_by_segment(BAF_heatmap, segments, outheatmap_dir_name, j)


def plot_BAF_heatmap_by_segment(BAF_heatmap, segments, outheatmap_dir_name, seg_idx):
    color_max = BAF_heatmap.color_max
    seg_name = segments[seg_idx][0]
    BAF_counts = BAF_heatmap.BAF_counts[seg_idx]
    
    pyplot.figure(figsize=(8,8))
    pyplot.xlim((0, 100))
    pyplot.ylim((0, 100))
    pyplot.xticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
    pyplot.yticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
    pyplot.xlabel('Tumour genome B allele frequency (BAF)')
    pyplot.ylabel('Normal genome B allele frequency (BAF)')
    pyplot.imshow(BAF_counts, vmin = 0, vmax = color_max)
    pyplot.colorbar(ticks=[0, color_max], orientation='vertical')
    pyplot.savefig('./' + outheatmap_dir_name + '/' + seg_name)
    