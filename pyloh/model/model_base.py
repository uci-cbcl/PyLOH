'''
Created on 2013-07-29

@author: Yi Li
'''
from ConfigParser import ConfigParser

import numpy as np

from pyloh import constants
from pyloh.preprocessing.io import Data

class ProbabilisticModel(object):
    def __init__(self):
        self.priors_parser = PriorParser()
        self.data = Data()
        self._init_components()
        
    
    def read_priors(self, priors_filename):
        self.priors_parser.read_priors(priors_filename)
        self.priors = self.priors_parser.priors
    
    def read_data(self, filename_base):
        self.data.read_data(filename_base)
        
    def preprocess_data(self):
        raise NotImplemented
        
    def run(self, max_iters, stop_value):
        trainer = self.model_trainer_class(self.priors, self.data, max_iters, stop_value)
        
        trainer.train()
        
        self.model_parameters = trainer.model_parameters
        
    def write_parameters(self, filename_base):
        self.model_parameters.write_parameters(filename_base)

    def _init_components(self):
        raise NotImplemented

#JointSNVMix
class ModelTrainer(object):
    def __init__(self, priors, data, max_iters, stop_value):
        self.priors = priors
        
        self.data = data
        
        self.max_iters = max_iters
        
        self.stop_value = stop_value
        
        self.iters = 0
            
        self._init_components()
        
    def train(self):
        converged = False
        
        parameters = self.model_parameters.parameters
        old_log_likelihood = self.log_likelihood.get_log_likelihood(parameters)
  
        while converged == False:
            self._E_step()
            self._M_step()

            log_likelihood = self.log_likelihood.get_log_likelihood(parameters)
            
            if self.iters > 0:
                ll_change = (log_likelihood - old_log_likelihood) / np.abs(old_log_likelihood)
            else:
                ll_change = float('inf')
            
            self._print_running_info(self.iters, log_likelihood, old_log_likelihood, ll_change)
            
            old_log_likelihood = log_likelihood
            
            if np.abs(ll_change) < self.stop_value:
                converged = True
                            
            if self.iters >= self.max_iters:
                print "Maximum numbers of EM iterations exceeded. Exiting training..."                
                converged = True
            
            self.iters += 1

    def _E_step(self):
        self.latent_variables.update(self.model_parameters.parameters, self.iters)
        
    def _M_step(self):
        self.model_parameters.update(self.latent_variables.sufficient_statistics)

    def _print_running_info(self, iters, log_likelihood, old_log_likelihood, ll_change):
        raise NotImplemented  

    def _init_components(self):
        raise NotImplemented

class LatentVariables(object): 
    def __init__(self, data):       
        self.data = data

    def update(self, parameters):
        raise NotImplemented

class ModelParameters(object):
    def __init__(self, priors, data):
        self.priors = priors
        self.data = data
        
        self._init_parameters()
    
    def update(self, sufficient_statistics):
        raise NotImplemented
    
    def write_parameters(filename_base):
        raise NotImplemented

    def _init_parameters(self):
        raise NotImplemented
    
class LogLikelihood(object):
    def __init__(self, data):
        self.data = data
        
    def get_log_likelihood(self, parameters):
        raise NotImplemented
    
#JointSNVMix        
class PriorParser(object):
    def __init__(self):                
        self.priors = {}
        
    def read_priors(self, priors_filename):                
        self.parser = ConfigParser()
        self.parser.read(priors_filename)
        
        self.priors['omega'] = np.zeros(constants.GENOTYPES_TUMOR_NUM)
        
        for i, genotype in enumerate(constants.GENOTYPES_TUMOR):
            self.priors['omega'][i] = self.parser.getfloat('omega', genotype)
    