'''
Created on 2013-07-31

@author: Yi Li
'''
import sys

import numpy as np

from pyloh import constants
from pyloh.preprocess.data import Data
from pyloh.model.model_base import *
from pyloh.model.utils import *

class PoissonProbabilisticModel(ProbabilisticModel):
    def __init__(self, allele_number_max):
        ProbabilisticModel.__init__(self, allele_number_max)

    def read_priors(self, priors_filename):
        if priors_filename != None:
            self.priors_parser.read_priors(priors_filename, self.allele_number_max)
            self.priors = self.priors_parser.priors
        else:
            self.priors = {}
            self.priors['omega'] = np.array(get_omega(self.allele_number_max))*1.0

    def preprocess(self):
        config_parameters = {}
        config_parameters['allele_number_max'] = self.allele_number_max
        config_parameters['genotypes_tumor'] =  get_genotypes_tumor(self.allele_number_max)
        config_parameters['genotypes_tumor_num'] = get_genotypes_tumor_num(self.allele_number_max)
        config_parameters['alleletypes_tumor'] = get_alleletypes_tumor(self.allele_number_max)
        config_parameters['copynumber_tumor'] = get_copynumber_tumor(self.allele_number_max)
        config_parameters['copynumber_tumor_num'] = get_copynumber_tumor_num(self.allele_number_max)
        config_parameters['MU_T'] = get_MU_T(self.allele_number_max)
        config_parameters['P_CG'] = get_P_CG(self.allele_number_max)
        
        self.config_parameters = config_parameters
        self.data.segments.compute_Lambda_S()
        
    def _init_components(self):
        self.model_trainer_class = PoissonModelTrainer
        
class PoissonModelTrainer(ModelTrainer):
    def __init__(self, priors, data, idx_restart, restart_parameters, config_parameters, max_iters, stop_value):
        ModelTrainer.__init__(self, priors, data, idx_restart, restart_parameters, config_parameters, max_iters, stop_value)
        
    def _init_components(self):
        self.latent_variables = PoissonLatentVariables(self.data, self.restart_parameters, self.config_parameters)
        
        self.model_parameters = PoissonModelParameters(self.priors, self.data, self.restart_parameters, self.config_parameters)
        
        self.model_likelihood = PoissonModelLikelihood(self.data, self.restart_parameters, self.config_parameters)
    
    def _print_running_info(self, new_log_likelihood, old_log_likelihood, ll_change):
        c_S = self.restart_parameters['copy_number_base']
        phi_init = self.restart_parameters['phi_init']
        allele_number_max = self.config_parameters['allele_number_max']
        
        print "#" * 100
        print "# Running Info."
        print "#" * 100
        print "Round of restarts : ", self.idx_restart + 1
        print "Baseline copy number : ", c_S
        print "Maximum copy number of each allele : ", allele_number_max
        print "Initial tumor celluarity : ", phi_init
        print "Number of iterations : ", self.iters
        print "New log-likelihood : ", new_log_likelihood
        print "Old log-likelihood : ", old_log_likelihood 
        print "Log-likelihood change : ", ll_change
        print "Parameters :"
        print "Tumor purity by CNV : {0:.3f}".format(self.model_parameters.parameters['phi_CNV'])
        print "Tumor purity by LOH : {0:.3f}".format(self.model_parameters.parameters['phi_LOH'])
        print "Tumor purity combined : {0:.3f}".format(self.model_parameters.parameters['phi'])
        sys.stdout.flush()
        
class PoissonLatentVariables(LatentVariables):
    def __init__(self, data, restart_parameters, config_parameters):
        LatentVariables.__init__(self, data, restart_parameters, config_parameters)
    
    def update(self, parameters, iters):
        eta = constants.ETA
        burn_in = constants.BURN_IN
        
        if iters >= burn_in:
            eta = 1
        
        J = self.data.seg_num
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        
        sufficient_statistics = {}
        
        sufficient_statistics['xi'] = []
        sufficient_statistics['psi'] =  np.zeros((J, C))
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                sufficient_statistics['xi'].append('NONE')
                continue
                        
            xi_j, psi_j = self._sufficient_statistics_by_segment(parameters, eta, j)
            sufficient_statistics['xi'].append(xi_j)
            sufficient_statistics['psi'][j, :] = psi_j
            
        self.sufficient_statistics = sufficient_statistics
    
    def _sufficient_statistics_by_segment(self, parameters, eta, j):        
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        I_j = d_T_j.shape[0]
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        rho_j = parameters['rho'][j]
        phi = parameters['phi']
        
        c_S = self.restart_parameters['copy_number_base']
        Lambda_S = self.data.segments.Lambda_S
        
        xi_j = self._get_xi_by_segment(b_T_j, d_T_j, I_j, phi, eta)
        
        kappa_j = self._get_kappa_by_segment(b_T_j, d_T_j, I_j, phi, xi_j, eta)
        
        psi_j = self._get_psi_by_segment(D_N_j, D_T_j, phi, rho_j, kappa_j, eta)
        
        #log_likelihoods_LOH = poisson_ll_func_LOH(b_T_j, d_T_j, rho_j, phi, self.config_parameters)
        #log_likelihoods_CNV = poisson_ll_func_CNV(D_N_j, D_T_j, c_S, Lambda_S, rho_j, phi, self.config_parameters)
        #log_likelihoods = poisson_ll_func(b_T_j, d_T_j, c_S, D_N_j, D_T_j, Lambda_S, rho_j, phi, config_parameters)

        #xi_j = log_space_normalise_rows_annealing(log_likelihoods_LOH, eta)
        #psi_j = log_space_normalise_rows_annealing(log_likelihoods_CNV, eta)
        #psi_j = log_space_normalise_rows_annealing(log_likelihoods, eta)
        
        #u_j = np.add.reduce(xi_j, axis = 0)
        #v_j = xi_j*np.dot(b_T_j.reshape(I_j, 1), np.ones((1, G)))
        #v_j = np.add.reduce(v_j, axis = 0)
        #w_j = xi_j*np.dot(d_T_j.reshape(I_j, 1), np.ones((1, G)))
        #w_j = np.add.reduce(w_j, axis = 0)
        #
        #return psi_j, u_j, v_j, w_j
        
        return xi_j, psi_j

    def _get_xi_by_segment(self, b_T_j, d_T_j, I_j, phi, eta):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        p_CG = self.config_parameters['P_CG']
        eps = constants.EPS
        
        xi_j = np.zeros((I_j, C, G))
        
        for c in range(0, C):
            mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_T[c], phi)
            log_likelihoods = np.log(p_CG[c]) + log_binomial_likelihood(b_T_j, d_T_j, mu_E_c)
            xi_j[:, c, :] = log_space_normalise_rows_annealing(log_likelihoods, eta)
        
        xi_j += eps
        
        return xi_j
    
    def _get_kappa_by_segment(self, b_T_j, d_T_j, I_j, phi, xi_j, eta):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        
        kappa_j = np.zeros(C)
        
        for c in range(0, C):
            mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_T[c], phi)
            log_likelihoods = np.log(xi_j[:, c, :]) + log_binomial_likelihood(b_T_j, d_T_j, mu_E_c)
            temp = np.logaddexp.reduce(log_likelihoods, axis = 1)
            kappa_j[c] = temp.sum()
        
        return kappa_j
    
    def _get_psi_by_segment(self, D_N_j, D_T_j, phi, rho_j, kappa_j, eta):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_S = self.restart_parameters['copy_number_base']
        Lambda_S = self.data.segments.Lambda_S
    
        c_E_j = get_c_E(c_N, c_T, phi)
        c_E_s = get_c_E(c_N, c_S, phi)
        lambda_E_j = D_N_j*c_E_j*Lambda_S/c_E_s
        log_likelihoods = np.log(rho_j) + log_poisson_likelihood(D_T_j, lambda_E_j) + kappa_j
        psi_j = log_space_normalise_rows_annealing(log_likelihoods, eta)
        
        return psi_j

class PoissonModelParameters(ModelParameters):
    def __init__(self, priors, data, restart_parameters, config_parameters):
        ModelParameters.__init__(self, priors, data, restart_parameters, config_parameters)
    
    def update(self, sufficient_statistics):
        self.sufficient_statistics = sufficient_statistics
        
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
        G = self.config_parameters['genotypes_tumor_num']
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
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])

        J = self.data.seg_num
        G = self.config_parameters['genotypes_tumor_num']
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
        
        phi_init = self.restart_parameters['phi_init']
        
        parameters['phi'] = phi_init
        parameters['phi_CNV'] = phi_init
        parameters['phi_LOH'] = phi_init
        parameters['rho'] = self._update_rho_by_priors()
        parameters['rho_CNV'] = parameters['rho']
        parameters['rho_LOH'] = parameters['rho']
        
        self.parameters = parameters
        
    def _update_rho_by_priors(self):
        J = self.data.seg_num
        G = self.config_parameters['genotypes_tumor_num']
        
        omega = self.priors['omega']
        
        rho = np.zeros((J, G))
        
        for j in range(0, J):
            rho[j] = omega/omega.sum()
            
        return rho
    
    def _update_phi_by_segment(self, D_N_j, D_T_j, mu_E_j, rho_j, rho_CNV_j, rho_LOH_j):
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        eps = constants.EPS
        
        G = self.config_parameters['genotypes_tumor_num']
        
        c_S = self.restart_parameters['copy_number_base']
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
    
    def write_parameters(self, filename_base):
        outpurity_file_name = filename_base + '.PyLOH.purity'
        outseg_ext_file_name = filename_base + '.PyLOH.segments.extended'
        
        self._write_purity(outpurity_file_name)
        self._write_seg_extended(outseg_ext_file_name)
        
    def _write_purity(self, outpurity_file_name):
        outfile = open(outpurity_file_name, 'w')
        
        c_S = self.restart_parameters['copy_number_base']
        allele_number_max = self.config_parameters['allele_number_max']
        
        outfile.write("Optimum baseline copy number : {0}".format(c_S) + '\n')
        outfile.write("Maximum copy number of each allele : {0}".format(allele_number_max) + '\n')
        outfile.write("Tumor purity by CNV : {0:.3f}".format(self.parameters['phi_CNV']) + '\n')
        outfile.write("Tumor purity by LOH : {0:.3f}".format(self.parameters['phi_LOH']) + '\n')
        outfile.write("Tumor purity combined : {0:.3f}".format(self.parameters['phi']) + '\n')
        
        outfile.close()
        
    def _write_seg_extended(self, outseg_ext_file_name):
        outfile = open(outseg_ext_file_name, 'w')

        outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num',
                                 'tumor_reads_num', 'LOH_frec', 'LOH_status', 'log2_ratio',
                                 'allele_type', 'copy_number']) + '\n')

        J = self.data.seg_num
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                allele_type_j = 'NONE'
                copy_number_j = 'NONE'
            else:
                allele_type_j = self._get_allele_type_by_segment(j)
                copy_number_j = self._get_copy_number_by_segment(j)
                
            segment_j = list(self.data.segments[j])
            segment_j.extend([allele_type_j, copy_number_j])
            
            outfile.write('\t'.join(map(str, segment_j)) + '\n')
        
        outfile.close()
        
    def _get_allele_type_by_segment(self, j):
        allele_types_tumor = self.config_parameters['alleletypes_tumor']
        G = self.config_parameters['genotypes_tumor_num']
        
        psi_j = self.sufficient_statistics['psi'][j]
        allele_types_prob = {}
        
        for g in range(0, G):
            allele_type = allele_types_tumor[g]
            if allele_type not in allele_types_prob.keys():
                allele_types_prob[allele_type] = psi_j[g]
            else:
                allele_types_prob[allele_type] += psi_j[g]
        
        allele_type_max = max(allele_types_prob, key=allele_types_prob.get)
        
        return allele_type_max

    def _get_copy_number_by_segment(self, j):
        copy_number_tumor = self.config_parameters['copynumber_tumor']
        G = self.config_parameters['genotypes_tumor_num']
        
        psi_j = self.sufficient_statistics['psi'][j]
        copy_number_prob = {}
        
        for g in range(0, G):
            copy_number = copy_number_tumor[g]
            if copy_number not in copy_number_prob.keys():
                copy_number_prob[copy_number] = psi_j[g]
            else:
                copy_number_prob[copy_number] += psi_j[g]
        
        copy_number_max = max(copy_number_prob, key=copy_number_prob.get)
        
        return copy_number_max

class PoissonModelLikelihood(ModelLikelihood):
    def __init__(self, data, restart_parameters, config_parameters):
        ModelLikelihood.__init__(self, data, restart_parameters, config_parameters)
    
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
        G = self.config_parameters['genotypes_tumor_num']
        
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        I_j = d_T_j.shape[0]
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        rho_j = parameters['rho'][j]
        phi = parameters['phi']
        
        Lambda_S = self.data.segments.Lambda_S
        c_S = self.restart_parameters['copy_number_base']
        
        log_likelihoods_LOH = poisson_ll_func_LOH(b_T_j, d_T_j, rho_j, phi, self.config_parameters)
        log_likelihoods_CNV = poisson_ll_func_CNV(D_N_j, D_T_j, c_S, Lambda_S, rho_j, phi, self.config_parameters)
        
        log_likelihoods_LOH = np.logaddexp.reduce(log_likelihoods_LOH, axis = 1)
        log_likelihoods_CNV = np.logaddexp.reduce(log_likelihoods_CNV, axis = 1)
        
        log_likelihood_j = log_likelihoods_LOH.sum() + log_likelihoods_CNV.sum()
        
        return log_likelihood_j
    

#===============================================================================
# Function
#===============================================================================
def poisson_ll_func(b_T, d_T, c_S, D_N, D_T, Lambda_S, rho, phi, config_parameters):
    I = d_T.shape[0]
    c_N = np.array(constants.COPY_NUMBER_NORMAL)
    c_T = np.array(config_parameters['copynumber_tumor'])
    mu_N = np.array(constants.MU_N)
    mu_T = np.array(config_parameters['MU_T'])
    c_E_j = get_c_E(c_N, c_T, phi)
    c_E_s = get_c_E(c_N, c_S, phi)
    
    lambda_E = D_N*c_E_j*Lambda_S/c_E_s
    
    mu_E = get_mu_E(mu_N, mu_T, c_N, c_T, phi)
        
    log_likelihoods = np.log(rho) + log_poisson_likelihood(D_T, lambda_E) \
    + np.add.reduce(log_binomial_likelihood(b_T, d_T, mu_E), axis = 0)
    
    return log_likelihoods

def poisson_ll_func_LOH(b_T, d_T, rho, phi, config_parameters):
    I = d_T.shape[0]
    c_N = np.array(constants.COPY_NUMBER_NORMAL)
    c_T = np.array(config_parameters['copynumber_tumor'])
    mu_N = np.array(constants.MU_N)
    mu_T = np.array(config_parameters['MU_T'])
    
    mu_E = get_mu_E(mu_N, mu_T, c_N, c_T, phi)
        
    log_likelihoods = np.log(rho) + log_binomial_likelihood(b_T, d_T, mu_E)
    
    return log_likelihoods
    
def poisson_ll_func_CNV(D_N, D_T, c_S, Lambda_S, rho, phi, config_parameters):
    c_N = np.array(constants.COPY_NUMBER_NORMAL)
    c_T = np.array(config_parameters['copynumber_tumor'])
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
            restart_parameters = {}
            restart_parameters['copy_number_base'] = c_S
            restart_parameters['phi_init'] = phi_init
            restart_parameters_list.append(restart_parameters)
    
    return restart_parameters_list




