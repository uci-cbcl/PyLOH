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
    def __init__(self, allelenumber_max):
        ProbabilisticModel.__init__(self, allelenumber_max)

    def read_priors(self, priors_filename):
        if priors_filename != None:
            self.priors_parser.read_priors(priors_filename, self.allelenumber_max)
            self.priors = self.priors_parser.priors
        else:
            self.priors = {}
            self.priors['omega'] = np.array(get_omega(self.allelenumber_max))*1.0

    def preprocess(self):
        config_parameters = {}
        config_parameters['allelenumber_max'] = self.allelenumber_max
        config_parameters['genotypes_tumor'] =  get_genotypes_tumor(self.allelenumber_max)
        config_parameters['genotypes_tumor_num'] = get_genotypes_tumor_num(self.allelenumber_max)
        config_parameters['copynumber_tumor'] = get_copynumber_tumor(self.allelenumber_max)
        config_parameters['copynumber_tumor_compat'] = get_copynumber_tumor_compat(self.allelenumber_max)
        config_parameters['copynumber_tumor_num'] = get_copynumber_tumor_num(self.allelenumber_max)
        config_parameters['MU_T'] = get_MU_T(self.allelenumber_max)
        config_parameters['P_CG'] = get_P_CG(self.allelenumber_max)
        
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
        allelenumber_max = self.config_parameters['allelenumber_max']
        
        print "#" * 100
        print "# Running Info."
        print "#" * 100
        print "Round of restarts : ", self.idx_restart + 1
        print "Baseline copy number : ", c_S
        print "Maximum copy number of each allele : ", allelenumber_max
        print "Initial tumor purity : ", phi_init
        print "Number of iterations : ", self.iters
        print "New log-likelihood : ", new_log_likelihood
        print "Old log-likelihood : ", old_log_likelihood 
        print "Log-likelihood change : ", ll_change
        print "Tumor purity : {0:.3f}".format(self.model_parameters.parameters['phi'])
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
        sufficient_statistics['kappa'] = np.zeros((J, C))
        sufficient_statistics['psi'] =  np.zeros((J, C))
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                sufficient_statistics['xi'].append('NONE')
                continue
                        
            xi_j, kappa_j, psi_j = self._sufficient_statistics_by_segment(parameters, eta, j)
            sufficient_statistics['xi'].append(xi_j)
            sufficient_statistics['kappa'][j, :] = kappa_j
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
        
        return xi_j, kappa_j, psi_j

    def _get_xi_by_segment(self, b_T_j, d_T_j, I_j, phi, eta):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_gT = np.array(self.config_parameters['copynumber_tumor_compat'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_gT, phi)
        p_CG = self.config_parameters['P_CG']
        eps = constants.EPS
        
        xi_j = np.zeros((I_j, C, G))
        
        for c in range(0, C):
            log_likelihoods = np.log(p_CG[c]) + log_binomial_likelihood(b_T_j, d_T_j, mu_E_c)
            xi_j[:, c, :] = log_space_normalise_rows_annealing(log_likelihoods, eta)
        
        xi_j += eps
        
        return xi_j
    
    def _get_kappa_by_segment(self, b_T_j, d_T_j, I_j, phi, xi_j, eta):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_gT = np.array(self.config_parameters['copynumber_tumor_compat'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_gT, phi)
        
        kappa_j = np.zeros(C)
        
        for c in range(0, C):
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
        
        xi = sufficient_statistics['xi']
        psi = sufficient_statistics['psi']
        
        rho = self._update_rho(psi)
        phi = self._update_phi(xi, psi)
        
        parameters['rho'] = rho
        parameters['phi'] = phi
        
        self.parameters = parameters
    
    def _update_rho(self, psi):
        x1 = constants.UPDATE_WEIGHTS['x1']
        y1 = constants.UPDATE_WEIGHTS['y1']

        rho_data = psi
        rho_priors = self._update_rho_by_priors()
        
        rho = x1*rho_data + y1*rho_priors
        
        return rho
    
    def _update_phi(self, xi, psi):
        phi_start = 0.01
        phi_end = 0.99
        phi_stop = 1e-5
        phi_change = 1
        
        while phi_change > phi_stop:
            phi_left = phi_start + (phi_end - phi_start)*1/3
            phi_right = phi_start + (phi_end - phi_start)*2/3
            
            unnorm_complete_ll_left = self._complete_log_likelihood(phi_left, xi, psi)
            unnorm_complete_ll_right = self._complete_log_likelihood(phi_right, xi, psi)
            
            if unnorm_complete_ll_left >= unnorm_complete_ll_right:
                phi_change = phi_end - phi_right
                phi_end = phi_right
            else:
                phi_change = phi_left - phi_start
                phi_start = phi_left
                            
        phi_optimum = (phi_start + phi_end)/2  
        
        return phi_optimum
    
    def _complete_log_likelihood(self, phi, xi, psi):
        J = self.data.seg_num
        
        complete_ll = 0
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                continue            

            complete_ll += self._complete_ll_CNA_by_segment(phi, psi[j], j)
            complete_ll += self._complete_ll_LOH_by_segment(phi, xi[j], psi[j], j)
        
        return complete_ll
    
    def _complete_ll_CNA_by_segment(self, phi, psi_j, j):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_S = self.restart_parameters['copy_number_base']
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        Lambda_S = self.data.segments.Lambda_S
    
        c_E_j = get_c_E(c_N, c_T, phi)
        c_E_s = get_c_E(c_N, c_S, phi)
        lambda_E_j = D_N_j*c_E_j*Lambda_S/c_E_s
        
        complete_ll_CNA_j = (psi_j*(D_T_j*np.log(lambda_E_j) - lambda_E_j)).sum()
        
        return complete_ll_CNA_j
    
    def _complete_ll_LOH_by_segment(self, phi, xi_j, psi_j, j):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_gT = np.array(self.config_parameters['copynumber_tumor_compat'])
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_gT, phi)
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        
        complete_ll_LOH_j = 0
        
        for c in range(0, C):
            log_likelihoods = log_binomial_likelihood(b_T_j, d_T_j, mu_E_c)
            log_likelihoods *= xi_j[:, c, :]
            
            complete_ll_LOH_j += psi_j[c]*log_likelihoods.sum()
        
        return complete_ll_LOH_j

    def _init_parameters(self):
        parameters = {}
        
        phi_init = self.restart_parameters['phi_init']
        
        parameters['phi'] = phi_init
        parameters['rho'] = self._update_rho_by_priors()
        
        self.parameters = parameters
        
    def _update_rho_by_priors(self):
        J = self.data.seg_num
        C = self.config_parameters['copynumber_tumor_num']
        
        omega = self.priors['omega']
        
        rho = np.zeros((J, C))
        
        for j in range(0, J):
            rho[j] = omega/omega.sum()
            
        return rho
       
    def write_parameters(self, filename_base):
        outpurity_file_name = filename_base + '.PyLOH.purity'
        outseg_ext_file_name = filename_base + '.PyLOH.segments.extended'
        
        self._write_purity(outpurity_file_name)
        self._write_seg_extended(outseg_ext_file_name)
        
    def _write_purity(self, outpurity_file_name):
        outfile = open(outpurity_file_name, 'w')
        
        c_S = self.restart_parameters['copy_number_base']
        allelenumber_max = self.config_parameters['allelenumber_max']
        
        outfile.write("Optimum baseline copy number : {0}".format(c_S) + '\n')
        outfile.write("Maximum copy number of each allele : {0}".format(allelenumber_max) + '\n')
        outfile.write("Tumor purity : {0:.3f}".format(self.parameters['phi']) + '\n')
        
        outfile.close()
        
    def _write_seg_extended(self, outseg_ext_file_name):
        outfile = open(outseg_ext_file_name, 'w')

        outfile.write('\t'.join(['#seg_name', 'chrom', 'start', 'end', 'normal_reads_num',
                                 'tumor_reads_num', 'LOH_frec', 'LOH_status', 'log2_ratio',
                                 'copy_number']) + '\n')

        J = self.data.seg_num
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                copy_number_j = 'NONE'
            else:
                copy_number_j = self._get_copy_number_by_segment(j)
                
            segment_j = list(self.data.segments[j])
            segment_j.extend([copy_number_j])
            
            outfile.write('\t'.join(map(str, segment_j)) + '\n')
        
        outfile.close()

    def _get_copy_number_by_segment(self, j):
        C = self.config_parameters['copynumber_tumor_num']
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        
        psi_j = self.sufficient_statistics['psi'][j]
        idx_c_T_max = psi_j.argmax()
        c_T_max = c_T[idx_c_T_max]
        
        return c_T_max

class PoissonModelLikelihood(ModelLikelihood):
    def __init__(self, data, restart_parameters, config_parameters):
        ModelLikelihood.__init__(self, data, restart_parameters, config_parameters)
    
    def get_log_likelihood(self, parameters, priors):
        J = self.data.seg_num
        
        log_likelihood = 0
        
        for j in range(0, J):
            LOH_status_j = self.data.segments[j][7]
            
            if LOH_status_j == 'NONE':
                continue
            
            log_likelihood_j = self._get_log_likelihood_by_segment(parameters, priors, j)
            log_likelihood += log_likelihood_j
        
        return log_likelihood
    
    def _get_log_likelihood_by_segment(self, parameters, priors, j):
        C = self.config_parameters['copynumber_tumor_num']
        G = self.config_parameters['genotypes_tumor_num']
        c_N = np.array(constants.COPY_NUMBER_NORMAL)
        c_T = np.array(self.config_parameters['copynumber_tumor'])
        c_gT = np.array(self.config_parameters['copynumber_tumor_compat'])
        c_S = self.restart_parameters['copy_number_base']
        a_T_j = self.data.paired_counts[j][:, 2]
        b_T_j = self.data.paired_counts[j][:, 3]
        d_T_j = a_T_j + b_T_j
        D_N_j = self.data.segments[j][4]
        D_T_j = self.data.segments[j][5]
        Lambda_S = self.data.segments.Lambda_S
        rho_j = parameters['rho'][j]
        phi = parameters['phi']
        omega = priors['omega']
        p_CG = self.config_parameters['P_CG']
        eps = constants.EPS
    
        c_E_j = get_c_E(c_N, c_T, phi)
        c_E_s = get_c_E(c_N, c_S, phi)
        lambda_E_j = D_N_j*c_E_j*Lambda_S/c_E_s
        
        mu_N = np.array(constants.MU_N)
        mu_T = np.array(self.config_parameters['MU_T'])
        mu_E_c = get_mu_E(mu_N, mu_T, c_N, c_gT, phi)

        #log-likelihood for LOH
        ll_LOH_j = np.zeros(C)
        
        for c in range(0, C):
            log_likelihoods = np.log(p_CG[c]) + log_binomial_likelihood(b_T_j, d_T_j, mu_E_c)
            log_likelihoods = np.logaddexp.reduce(log_likelihoods, axis = 1)            
            ll_LOH_j[c] = log_likelihoods.sum()
            
        #log-likelihood for CNA
        ll_CNA_j = log_poisson_likelihood(D_T_j, lambda_E_j)
        
        ll_j = np.log(rho_j) + ll_LOH_j + ll_CNA_j
        
        ll_j = np.logaddexp.reduce(ll_j, axis = 1)[0]
        
        #Dirichlet priors for each segment
        ll_j = ll_j + log_dirichlet_pdf(rho_j, omega)
        
        return ll_j
    

#===============================================================================
# Function
#===============================================================================
   
def poisson_restart_parameters_list():
    restart_parameters_list = []
    
    for c_S in constants.COPY_NUMBER_BASE:
        for phi_init in constants.PHI_INIT:
            restart_parameters = {}
            restart_parameters['copy_number_base'] = c_S
            restart_parameters['phi_init'] = phi_init
            restart_parameters_list.append(restart_parameters)
    
    return restart_parameters_list




