'''
Created on 2012-08-13

@author: Yi Li
'''
import sys

import numpy as np

from pyloh import constants
from pyloh.model.poisson_model import PoissonProbabilisticModel, poisson_restart_parameters_list

def run_poisson_model(args):
    model_parameters_list = []
    log_likelihood_list = []
    restart_parameters_list = poisson_restart_parameters_list()
    restart_num = len(restart_parameters_list)
    
    model = PoissonProbabilisticModel(args.allele_number_max)
    model.read_priors(args.priors_file_name)
    model.read_data(args.filename_base)
    model.preprocess()

    for idx_restart in range(0, restart_num):
        restart_parameters = restart_parameters_list[idx_restart]
        model.run(idx_restart, restart_parameters, args.max_iters, args.stop_value)
        
        model_parameters_list.append(model.model_parameters)
        log_likelihood_list.append(model.log_likelihood)
        
    log_likelihood_list = np.array(log_likelihood_list)
    
    idx_restart_optimum = log_likelihood_list.argmax()
    
    model_parameters_optimum = model_parameters_list[idx_restart_optimum]
    restart_parameters_optimum = restart_parameters_list[idx_restart_optimum]
    c_S_optimum = restart_parameters_optimum['copy_number_base']
    model_parameters_optimum.write_parameters(args.filename_base)
    
    print "*" * 100
    print "* Finish."
    print "*" * 100
    print "Optimum log-likelihood : ", log_likelihood_list[idx_restart_optimum]
    print "Optimum baseline copy number : ", c_S_optimum
    print "Maximum copy number of each allele : ", args.allele_number_max
    print "Tumor cellular frequency by CNV : {0:.3f}".format(model_parameters_optimum.parameters['phi_CNV'])
    print "Tumor cellular frequency by LOH : {0:.3f}".format(model_parameters_optimum.parameters['phi_LOH'])
    print "Tumor cellular frequency combined : {0:.3f}".format(model_parameters_optimum.parameters['phi'])
    sys.stdout.flush()
    
