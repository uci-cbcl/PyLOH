'''
Created on 2013-07-31

@author: Yi Li
'''
import numpy as np

from pyloh import constants
from pyloh.preprocessing.io import Data
from pyloh.model.model_base import *
from pyloh.model.utils import *

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
        
        self.model_parameters = PoissonModelParameters(self.priors, self.data)
        
        self.log_likelihood = PoissonLogLikelihood(self.data)
    
    def _print_running_info(self, iters, log_likelihood, old_log_likelihood, ll_change):
        print "#" * 100
        print "# Running Info."
        print "#" * 100
        print "Number of iterations : ", iters
        print "New log-likelihood : ", log_likelihood
        print "Old log-likelihood : ", old_log_likelihood 
        print "Log-likelihood change : ", ll_change
    
        print "Parameters :"
        print "Tumor celluar frequency: {0}".format(self.model_parameters.parameters['phi']) 
        
class PoissonLatentVariables(LatentVariables):
    def __init__(self, data):
        LatentVariables.__init__(self, data)
    
    def update(self, parameters, iters):
        eta = constants.ETA
        burn_in = constants.BURN_IN
        
        if iters >= burn_in:
            eta = 1
        
        J = self.data.seg_num
        G = constants.GENOTYPES_TUMOR_NUM
        
        sufficient_statistics = {}
        
        sufficient_statistics['psi'] = np.zeros((J, G))
        sufficient_statistics['u'] = np.zeros((J, G))
        sufficient_statistics['v'] = np.zeros((J, G))
        sufficient_statistics['w'] = np.zeros((J, G))
        
        for j in range(0, J):
            psi_j, u_j, v_j, w_j = self._sufficient_statistics_by_segment(parameters, eta, j)
            sufficient_statistics['psi'][j, :] = psi_j
            sufficient_statistics['u'][j, :] = u_j
            sufficient_statistics['v'][j, :] = v_j
            sufficient_statistics['w'][j, :] = w_j
            
        self.sufficient_statistics = sufficient_statistics
    
    def _sufficient_statistics_by_segment(self, parameters, eta, j):        
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        I_j = d_T_j.shape[0]
        Lambda_S = self.data.segments.Lambda_S
        
        log_likelihoods_paired_counts = poisson_ll_func_paired_counts(a_T_j, b_T_j, d_T_j, parameters)
        log_likelihoods_reads_depth = poisson_ll_func_reads_depth(D_N_j, D_T_j, Lambda_S, parameters)

        xi_j = log_space_normalise_rows_annealing(log_likelihoods_paired_counts, eta)
        psi_j = log_space_normalise_rows_annealing(log_likelihoods_reads_depth, eta)
        
        u_j = np.add.reduce(xi_j, axis = 0)
        v_j = xi_j*np.dot(b_T_j.reshape(I_j, 1), np.ones((1, G)))
        v_j = np.add.reduce(v_j, axis = 0)
        w_j = xi_j*np.dot(d_T_j.reshape(I_j, 1), np.ones((1, G)))
        w_j = np.add.reduce(w_j, axis = 0)
        
        return psi_j, u_j, v_j, w_j

#===============================================================================
# Function
#===============================================================================
def poisson_ll_func_paired_counts(a_T, b_T, d_T, parameters):
    
def poisson_ll_func_reads_depth(D_N, D_T, Lambda_S, parameters):




