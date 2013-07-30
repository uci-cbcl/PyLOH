'''
Created on 2013-07-29

@author: Yi Li
'''
from ConfigParser import ConfigParser

import numpy as np

from pyloh import constants
from pyloh.preprocessing.io import Data

class ProbabilisticModel(object):
    def __init__(self, model_trainer_class):
        self.priors_parser = PriorParser()
        self.model_trainer_class = model_trainer_class
        self.data = Data()
    
    def read_priors(self, priors_filename):
        self.priors_parser.read_priors(priors_filename)
        self.priors = self.priors_parser.priors
    
    def read_data(self, filename_base):
        self.data.read_data(filename_base)
        
    def run(self, max_iters, stop_value):
        trainer = self.model_trainer_class(self.priors, self.data, max_iters, stop_value)
        
        trainer.train(priors, data, max_iters, stop_value)
        
        self.parameters = trainer.parameters
        
    def write_parameters(self, filename_base):
        self.parameters.write_parameters(filename_base)

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

class ModelTrainer(object):
    def __init__(self, priors, data, max_iters, stop_value):



