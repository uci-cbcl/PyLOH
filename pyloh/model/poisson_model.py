'''
Created on 2013-07-31

@author: Yi Li
'''
import numpy as np

from pyloh import constants
from pyloh.preprocessing.io import Data
from pyloh.model.model_base import *

class PoissonProbabilisticModel(ProbabilisticModel):
    def __init__(self):
        ProbabilisticModel.__init__(self)
        
    def preprocess_data(self):
        self.data.segments.compute_Lambda_S()
        
    def _init_components(self):
        self.trainer_class = PoissonModelTrainer
        
class PoissonModelTrainer(ModelTrainer):
    def __init__(self, priors, data, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, max_iters, stop_value)
        
    def _init_components(self):
        self.latent_variables = PoissonLatentVariables(self.data)
        
        self.model_parameters = PoissonModelParameters(self.priors, self.data,
                                self.latent_variables.sufficient_statistics)
        
        self.log_likelihood = PoissonLogLikelihood(self.data)
        
class PoissonLatentVariables(LatentVariables):
    def __init__(self, data):
        LatentVariables.__init__(self, data)
    
    #def update(self, parameters):
    
    def _init_sufficient_statistics(self, data):
        sufficient_statistics = {}
        
        J = data.seg_num
        G = constants.GENOTYPES_TUMOR_NUM
        I = data.sites_num
        
        
        
