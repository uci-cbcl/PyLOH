'''
Created on 2012-08-18

@author: Yi Li
'''
import os
import sys

import numpy as np
import scipy

PLT_AVAIL = True
try:
    from matplotlib import pyplot
except:
    PLT_AVAIL = False

from pyloh import constants
from pyloh.preprocess.data import Segments, BAFHeatMap

def plot_BAF_heatmap(args):
    if PLT_AVAIL == True:
        pass
    else:
        print "matplotlib.pyplot not available, skip plotting..."
        sys.stdout.flush()
        sys.exit(-1)

    BAF_heatmap = BAFHeatMap()
    BAF_heatmap.read_heatmap(args.filename_base)
    BAF_heatmap.get_color_max()
    
    segments = Segments()
    inseg_file_name = args.filename_base + '.PyLOH.segments'
    segments.read_segfile(inseg_file_name)
    seg_num = segments.num
    
    outheatmap_dir_name = args.filename_base + '.PyLOH.heatmap.plot'
    if os.path.exists(outheatmap_dir_name) == False:
        os.mkdir(outheatmap_dir_name)
    
    for j in range(0, seg_num):
        BAF_counts_j = BAF_heatmap.BAF_counts[j]
        seg_name_j = segments[j][0]
        
        print 'Plotting segment {0}...'.format(seg_name_j)
        sys.stdout.flush()
        
        plot_BAF_heatmap_by_segment(BAF_counts_j, seg_name_j, outheatmap_dir_name)


def plot_BAF_heatmap_by_segment(BAF_counts, seg_name, outheatmap_dir_name):
    BAF_counts_min = constants.BAF_COUNTS_MIN
    BAF_counts_max = constants.BAF_COUNTS_MAX
    
    BAF_counts_sub = BAF_counts[BAF_counts_min:BAF_counts_max, BAF_counts_min:BAF_counts_max]
    color_max = BAF_counts_sub.max()
    
    pyplot.figure(figsize=(8,8), dpi = 150)
    pyplot.xlim((0, 100))
    pyplot.ylim((0, 100))
    pyplot.xticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
    pyplot.yticks(scipy.linspace(0, 100, 11), scipy.linspace(0, 1, 11))
    pyplot.xlabel('Tumor sample B allele frequency')
    pyplot.ylabel('Normal sample B allele frequency')
    pyplot.imshow(BAF_counts, vmin = 0, vmax = max(1, color_max))
    cbar = pyplot.colorbar(ticks=[0, color_max], orientation='vertical', shrink=0.78)
    cbar.ax.set_yticklabels(['0', '>= ' + str(int(color_max))])
    pyplot.savefig('./' + outheatmap_dir_name + '/' + seg_name, bbox_inches='tight')
    
