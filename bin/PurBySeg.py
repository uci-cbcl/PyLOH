#!/usr/bin/env python

import sys
import numpy as np
from pyloh import constants
from pyloh.preprocess.data import Data, Segments
from pyloh.model.utils import *
from pyloh.model.model_base import *
from pyloh.model.poisson_model import *

constants.COPY_NUMBER_BASE = [2]

def get_data_by_seg(data, j):
    
    segments = data.segments
    segments_j = Segments()
    segments_j.num = 1
    segments_j.Lambda_S = segments.Lambda_S
    segments_j.names.append(segments.names[j])
    segments_j.chroms.append(segments.chroms[j])
    segments_j.starts.append(segments.starts[j])
    segments_j.ends.append(segments.ends[j])
    segments_j.normal_reads_num.append(segments.normal_reads_num[j])
    segments_j.tumor_reads_num.append(segments.tumor_reads_num[j])
    segments_j.LOH_frec.append(segments.LOH_frec[j])
    segments_j.LOH_status.append(segments.LOH_status[j])
    segments_j.log2_ratio.append(segments.log2_ratio[j])
    
    data_j = Data()
    data_j.segments = segments_j
    data_j.paired_counts.append(data.paired_counts[j])
    data_j.sites_num.append(data.sites_num[j])
    data_j.seg_num = 1
    
    return data_j
    

def est_by_seg(data, config_parameters, j):
    priors = {}
    priors['omega'] = np.array(get_omega(config_parameters['allelenumber_max']))*1.0
    data_j = get_data_by_seg(data, j)
    model = PoissonProbabilisticModel(config_parameters['allelenumber_max'])
    model.priors = priors
    model.data = data_j
    model.config_parameters = config_parameters
    
    restart_num = len(constants.PHI_INIT)
    model_parameters_list = []
    log_likelihood_list = []
    for idx_restart in range(0, restart_num):
        phi_init = constants.PHI_INIT[idx_restart]
        
        restart_parameters = {}
        restart_parameters['copy_number_base'] = 2
        restart_parameters['phi_init'] = phi_init
        
        model.run(idx_restart, restart_parameters, 100, 1e-7)
        
        model_parameters_list.append(model.model_parameters)
        log_likelihood_list.append(model.log_likelihood)

    log_likelihood_list = np.array(log_likelihood_list)
    idx_restart_optimum = log_likelihood_list.argmax()
    model_parameters_optimum = model_parameters_list[idx_restart_optimum]
    purity_j = model_parameters_optimum.parameters['phi']
    copy_number_j = model_parameters_optimum._get_copy_number_by_segment(0)
    
    return purity_j, copy_number_j

def main():
    filename_base = sys.argv[1]
    outfile_name = filename_base + '.PyLOH.subclonal'
    outfile = open(outfile_name, 'w')
    
    data = Data()
    data.read_data(filename_base)
    data.segments.compute_Lambda_S()
    segments = data.segments
    
    allelenumber_max = 2 
    config_parameters = {}
    config_parameters['allelenumber_max'] = allelenumber_max
    config_parameters['genotypes_tumor'] =  get_genotypes_tumor(allelenumber_max)
    config_parameters['genotypes_tumor_num'] = get_genotypes_tumor_num(allelenumber_max)
    config_parameters['copynumber_tumor'] = get_copynumber_tumor(allelenumber_max)
    config_parameters['copynumber_tumor_compat'] = get_copynumber_tumor_compat(allelenumber_max)
    config_parameters['copynumber_tumor_num'] = get_copynumber_tumor_num(allelenumber_max)
    config_parameters['MU_T'] = get_MU_T(allelenumber_max)
    config_parameters['P_CG'] = get_P_CG(allelenumber_max)
    
    J = data.seg_num
    C = config_parameters['copynumber_tumor_num']
    G = config_parameters['genotypes_tumor_num']
    purity = np.ones(J)*-1
    copy_number = np.ones(J)*-1
    
    for j in range(0, J):
        LOH_status_j = data.segments[j][7]
        
        if LOH_status_j == 'NONE':
            continue
        
        purity_j, copy_number_j = est_by_seg(data, config_parameters, j)
        
        purity[j] = purity_j
        copy_number[j] = copy_number_j
        
    
    print '\t'.join(['#segment', 'chrom', 'LOH_frec', 'LOH_status', 'log2_ratio', 'purity', 'copy_number'])
    outfile.write('\t'.join(['#segment', 'chrom', 'LOH_frec', 'LOH_status', 'log2_ratio', 'purity', 'copy_number']) + '\n')
    for j in range(0, J):
        print '\t'.join(map(str, [segments[j][0], segments[j][1], segments[j][6], segments[j][7], segments[j][8], purity[j], int(copy_number[j])]))
        outfile.write('\t'.join(map(str, [segments[j][0], segments[j][1], segments[j][6], segments[j][7], segments[j][8], purity[j], int(copy_number[j])])) + '\n')
    



if __name__ == '__main__':
    main()