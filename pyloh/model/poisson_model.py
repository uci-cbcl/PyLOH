'''
Created on 2013-07-31

@author: Yi Li
'''
import numpy as np

from pyloh import constants
from pyloh.preprocess.data import Data
from pyloh.model.model_base import *
from pyloh.model.utils import *

class PoissonProbabilisticModel(ProbabilisticModel):
    def __init__(self):
        ProbabilisticModel.__init__(self)

    def read_priors(self, priors_filename):
        if priors_filename != None:
            self.priors_parser.read_priors(priors_filename)
            self.priors = self.priors_parser.priors
        else:
            self.priors = {}
            self.priors['omega'] = np.array(constants.OMEGA)*1.0

    def preprocess_data(self):
        self.data.segments.compute_Lambda_S()
        
    def _init_components(self):
        self.model_trainer_class = PoissonModelTrainer
        
class PoissonModelTrainer(ModelTrainer):
    def __init__(self, priors, data, idx_restart, restart_parameters, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, idx_restart, restart_parameters, max_iters, stop_value)
        
    def _init_components(self):
        self.latent_variables = PoissonLatentVariables(self.data, self.restart_parameters)
        
        self.model_parameters = PoissonModelParameters(self.priors, self.data, self.restart_parameters)
        
        self.model_likelihood = PoissonModelLikelihood(self.data, self.restart_parameters)
    
    def _print_running_info(self, idx_restart, restart_parameters, iters, new_log_likelihood, old_log_likelihood, ll_change):
        c_S, phi_init = restart_parameters
        
        print "#" * 100
        print "# Running Info."
        print "#" * 100
        print "Round of restarts : ", idx_restart + 1
        print "Baseline copy number : ", c_S
        print "Initial tumor celluar frequency : ", phi_init
        print "Number of iterations : ", iters
        print "New log-likelihood : ", new_log_likelihood
        print "Old log-likelihood : ", old_log_likelihood 
        print "Log-likelihood change : ", ll_change
        print "Parameters :"
        print "Tumor celluar frequency by CNV: {0}".format(self.model_parameters.parameters['phi_CNV'])
        print "Tumor celluar frequency by LOH: {0}".format(self.model_parameters.parameters['phi_LOH'])
        print "Tumor celluar frequency combined: {0}".format(self.model_parameters.parameters['phi']) 
        
class PoissonLatentVariables(LatentVariables):
    def __init__(self, data, restart_parameters):
        LatentVariables.__init__(self, data, restart_parameters)
    
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
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                continue
            
            psi_j, u_j, v_j, w_j = self._sufficient_statistics_by_segment(parameters, eta, j)
            sufficient_statistics['psi'][j, :] = psi_j
            sufficient_statistics['u'][j, :] = u_j
            sufficient_statistics['v'][j, :] = v_j
            sufficient_statistics['w'][j, :] = w_j
            
        self.sufficient_statistics = sufficient_statistics
    
    def _sufficient_statistics_by_segment(self, parameters, eta, j):
        G = constants.GENOTYPES_TUMOR_NUM
        
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        I_j = d_T_j.shape[0]
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        rho_j = parameters['rho'][j]
        phi = parameters['phi']
        
        c_S, __ = self.restart_parameters
        Lambda_S = self.data.segments.Lambda_S
        
        log_likelihoods_LOH = poisson_ll_func_LOH(b_T_j, d_T_j, rho_j, phi)
        log_likelihoods_CNV = poisson_ll_func_CNV(D_N_j, D_T_j, c_S, Lambda_S, rho_j, phi)

        xi_j = log_space_normalise_rows_annealing(log_likelihoods_LOH, eta)
        psi_j = log_space_normalise_rows_annealing(log_likelihoods_CNV, eta)
        
        u_j = np.add.reduce(xi_j, axis = 0)
        v_j = xi_j*np.dot(b_T_j.reshape(I_j, 1), np.ones((1, G)))
        v_j = np.add.reduce(v_j, axis = 0)
        w_j = xi_j*np.dot(d_T_j.reshape(I_j, 1), np.ones((1, G)))
        w_j = np.add.reduce(w_j, axis = 0)
        
        return psi_j, u_j, v_j, w_j

class PoissonModelParameters(ModelParameters):
    def __init__(self, priors, data, restart_parameters):
        ModelParameters.__init__(self, priors, data, restart_parameters)
    
    def update(self, sufficient_statistics):
        parameters = {}
                
        psi = sufficient_statistics['psi']
        u = sufficient_statistics['u']
        v = sufficient_statistics['v']
        w = sufficient_statistics['w']
        
        rho, rho_CNV, rho_LOH = self._update_rho(psi, u)
        phi, phi_CNV, phi_LOH = self._update_phi(v, w, rho, rho_CNV, rho_LOH)
        
        parameters['rho'] = rho
        parameters['rho_CNV'] = rho_CNV
        parameters['rho_LOH'] = rho_LOH
        parameters['phi'] = phi
        parameters['phi_CNV'] = phi_CNV
        parameters['phi_LOH'] = phi_LOH
        
        self.parameters = parameters
    
    def _update_rho(self, psi, u):
        x1 = constants.UPDATE_WEIGHTS['x1']
        y1 = constants.UPDATE_WEIGHTS['y1']
        z1 = constants.UPDATE_WEIGHTS['z1']
        
        J = self.data.seg_num
        G = constants.GENOTYPES_TUMOR_NUM
        I = np.array(self.data.sites_num)
        eps = constants.EPS
        
        rho_CNV = psi
        
        rho_LOH = u/np.dot((I + eps).reshape(J, 1), np.ones((1, G)))
        
        rho_priors = self._update_rho_by_priors()
        
        rho = x1*rho_CNV + y1*rho_LOH + z1*rho_priors
        
        return rho, rho_CNV, rho_LOH
    
    def _update_phi(self, v, w, rho, rho_CNV, rho_LOH):
        x2 = constants.UPDATE_WEIGHTS['x2']
        y2 = constants.UPDATE_WEIGHTS['y2']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(constants.COPY_NUMBER_TUMOR)
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(constants.MU_T)

        J = self.data.seg_num
        G = constants.GENOTYPES_TUMOR_NUM
        eps = constants.EPS

        phi_CNV = np.zeros(J)
        phi_LOH = np.zeros(J)
        weights_CNV = np.zeros(J)
        weights_LOH = np.zeros(J)
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                continue

            D_N_j = self.data.segments[j][4]
            D_T_j = self.data.segments[j][5]
            mu_E_j = v[j]/(w[j] + eps)
            
            for g in range(0, G):
                if mu_E_j[g] < min(mu_T):
                    mu_E_j[g] = min(mu_T)
                if mu_E_j[g] > max(mu_T):
                    mu_E_j[g] = max(mu_T)
            
            phi_CNV[j], phi_LOH[j], weights_CNV[j], weights_LOH[j] = self._update_phi_by_segment(D_N_j, D_T_j, mu_E_j,
                                                                                rho[j], rho_CNV[j], rho_LOH[j])
        
        phi_CNV_mean = np.average(phi_CNV, weights = weights_CNV)
        phi_LOH_mean = np.average(phi_LOH, weights = weights_LOH)
        
        phi = x2*phi_CNV_mean + y2*phi_LOH_mean
        
        return phi, phi_CNV_mean, phi_LOH_mean

    def _init_parameters(self):
        parameters = {}
        
        __, phi_init = self.restart_parameters
        
        parameters['phi'] = phi_init
        parameters['phi_CNV'] = phi_init
        parameters['phi_LOH'] = phi_init
        parameters['rho'] = self._update_rho_by_priors()
        parameters['rho_CNV'] = parameters['rho']
        parameters['rho_LOH'] = parameters['rho']
        
        self.parameters = parameters
        
    def _update_rho_by_priors(self):
        J = self.data.seg_num
        G = constants.GENOTYPES_TUMOR_NUM
        
        omega = self.priors['omega']
        
        rho = np.zeros((J, G))
        
        for j in range(0, J):
            rho[j] = omega/omega.sum()
            
        return rho
    
    def _update_phi_by_segment(self, D_N_j, D_T_j, mu_E_j, rho_j, rho_CNV_j, rho_LOH_j):
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(constants.COPY_NUMBER_TUMOR)
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(constants.MU_T)
        eps = constants.EPS
        
        G = constants.GENOTYPES_TUMOR_NUM
        
        c_S, __ = self.restart_parameters
        Lambda_S = self.data.segments.Lambda_S
        
        phi_CNV_j = 0
        phi_LOH_j = 0
        prob_sum_CNV_j = 0
        prob_sum_LOH_j = 0
         
        for g in range(0, G):
            phi_CNV_j_g = -1
            phi_LOH_j_g = -1
            
            if c_N[0] != c_T[g]:
                phi_CNV_j_g = get_phi_CNV(c_N[0], c_T[g], D_N_j, D_T_j, c_S, Lambda_S)
            
            if mu_N[0] != mu_T[g]:
                phi_LOH_j_g = get_phi_LOH(mu_N[0], mu_T[g], mu_E_j[g], c_N[0], c_T[g])
                
            if phi_CNV_j_g >= 0 and phi_CNV_j_g <= 1:
                phi_CNV_j = phi_CNV_j + rho_CNV_j[g]*phi_CNV_j_g
                prob_sum_CNV_j = prob_sum_CNV_j + rho_CNV_j[g]
                
            if phi_LOH_j_g >= 0 and phi_LOH_j_g <= 1:
                phi_LOH_j = phi_LOH_j + rho_LOH_j[g]*phi_LOH_j_g
                prob_sum_LOH_j = prob_sum_LOH_j + rho_LOH_j[g]
                
        if prob_sum_CNV_j == 0:
            prob_sum_CNV_j = eps
        if prob_sum_LOH_j == 0:
            prob_sum_LOH_j = eps
        
        phi_CNV_j = phi_CNV_j/prob_sum_CNV_j
        phi_LOH_j = phi_LOH_j/prob_sum_LOH_j
                
        return (phi_CNV_j, phi_LOH_j, prob_sum_CNV_j, prob_sum_LOH_j)

class PoissonModelLikelihood(ModelLikelihood):
    def __init__(self, data, restart_parameters):
        ModelLikelihood.__init__(self, data, restart_parameters)
    
    def get_log_likelihood(self, parameters):
        J = self.data.seg_num
        
        log_likelihood = 0
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                continue
            
            log_likelihood_j = self._get_log_likelihood_by_segment(parameters, j)
            log_likelihood += log_likelihood_j
            
        return log_likelihood
    
    def _get_log_likelihood_by_segment(self, parameters, j):
        G = constants.GENOTYPES_TUMOR_NUM
        
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        I_j = d_T_j.shape[0]
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        rho_j = parameters['rho'][j]
        phi = parameters['phi']
        
        Lambda_S = self.data.segments.Lambda_S
        c_S, __ = self.restart_parameters
        
        log_likelihoods_LOH = poisson_ll_func_LOH(b_T_j, d_T_j, rho_j, phi)
        log_likelihoods_CNV = poisson_ll_func_CNV(D_N_j, D_T_j, c_S, Lambda_S, rho_j, phi)
        
        log_likelihoods_LOH = np.logaddexp.reduce(log_likelihoods_LOH, axis = 1)
        log_likelihoods_CNV = np.logaddexp.reduce(log_likelihoods_CNV, axis = 1)
        
        log_likelihood_j = log_likelihoods_LOH.sum() + log_likelihoods_CNV.sum()
        
        return log_likelihood_j
    

#===============================================================================
# Function
#===============================================================================

def poisson_ll_func_LOH(b_T, d_T, rho, phi):
    I = d_T.shape[0]
    mu_N = np.array(constants.MU_N)
    mu_T = np.array(constants.MU_T)
    c_N = np.array(constants.COPY_NUMBER_NORMAL)
    c_T = np.array(constants.COPY_NUMBER_TUMOR)
    
    mu_E = get_mu_E(mu_N, mu_T, c_N, c_T, phi)
        
    log_likelihoods = np.log(rho) + log_binomial_likelihood(b_T, d_T, mu_E)
    
    return log_likelihoods
    
def poisson_ll_func_CNV(D_N, D_T, c_S, Lambda_S, rho, phi):
    c_N = np.array(constants.COPY_NUMBER_NORMAL)
    c_T = np.array(constants.COPY_NUMBER_TUMOR)
    c_S = np.array(c_S)
    
    c_E_j = get_c_E(c_N, c_T, phi)
    c_E_S = get_c_E(c_N, c_S, phi)
    
    lambda_E = D_N*c_E_j*Lambda_S/c_E_S
    
    log_likelihoods = np.log(rho) + log_poisson_likelihood(D_T, lambda_E)
    
    return log_likelihoods
    
def poisson_restart_parameters_list():
    restart_parameters_list = []
    
    for c_S in constants.COPY_NUMBER_BASE:
        for phi_init in constants.PHI_INIT:
            restart_parameters = (c_S, phi_init)
            restart_parameters_list.append(restart_parameters)
    
    return restart_parameters_list




