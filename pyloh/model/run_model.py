'''
Created on 2012-08-13

@author: Yi Li
'''
import numpy as np

from pyloh import constants
from pyloh.model.poisson_model import PoissonProbabilisticModel

def run_model(args):
    restart_num = constants.RESTART_NUM
    
    model_parameters_list = []
    log_likelihood_list = []
    
    model = PoissonProbabilisticModel()
    model.read_priors(args.priors_file_name)
    model.read_data(args.data_file_basename)
    model.preprocess_data()
    
    for i in range(0, restart_num):
        model.run(i, args.max_iters, args.stop_value)
        
        model_parameters_list.append(model.model_parameters)
        log_likelihood_list.append(model.log_likelihood)
        
    log_likelihood_list = np.array(log_likelihood_list)
    
    idx_restart_optimum = log_likelihood_list.argmax()
    
    model_parameters_optimum = model_parameters_list[idx_restart_optimum]
    
    print "*" * 100
    print "* Finish."
    print "*" * 100
    print "Optimum log-likelihood : ", log_likelihood_list[idx_restart_optimum]
    print "Tumor celluar frequency by CNV: {0}".format(model_parameters_optimum.parameters['phi_CNV'])
    print "Tumor celluar frequency by LOH: {0}".format(model_parameters_optimum.parameters['phi_LOH'])
    print "Tumor celluar frequency combined: {0}".format(model_parameters_optimum.parameters['phi']) 
    
