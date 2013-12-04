#!/usr/bin/env python

import numpy as np
import scipy
import sys
import os
from matplotlib import pyplot
from joint_snv_mix.file_formats.jcnt import JointCountsReader
from joint_snv_mix.classification.utils.log_pdf import log_binomial_likelihood
#from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows_annealing
from joint_snv_mix.classification.utils.normalise import log_space_normalise_rows

COUNTS_MIN = 10
COUNTS_MAX = 95
BAF_N_MIN = 0.4
BAF_N_MAX = 0.6
BAF_T_MIN = 0.35
BAF_T_MAX = 0.65
LOH_FREC_THRED = 0.25
SITES_NUM_MIN = 500
ERR = 0.01
EMPIRI_BAF = 0.485
EMPIRI_AAF = 1 - EMPIRI_BAF
PHI_INIT = 0.1
#ETA = 1.01
ETA = 1
BURN_IN = 30

GENOTYPES_NORMAL = ['AB']
GENOTYPES_TUMOR = ['A', 'B', 'AA', 'AB', 'BB', 'AAB', 'ABB', 'AABB']
GENOTYPES_TUMOR_NUM = len(GENOTYPES_TUMOR)
C_N = [2]
C_T = [1, 1, 2, 2, 2, 3, 3, 4]
MU_N = [EMPIRI_BAF/(EMPIRI_BAF + EMPIRI_AAF)]
MU_T = [ERR, 1-ERR, ERR, EMPIRI_BAF/(EMPIRI_BAF + EMPIRI_AAF), 1-ERR,
        EMPIRI_BAF*1/(EMPIRI_BAF*1 + EMPIRI_AAF*2),
        EMPIRI_BAF*2/(EMPIRI_BAF*2 + EMPIRI_AAF*1),
        EMPIRI_BAF/(EMPIRI_BAF + EMPIRI_AAF)]
omega = np.array([1000, 1000, 2, 100, 2, 2, 2, 2])


CHROM_LIST = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
              'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
              'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

def normal_heterozygous_filter(counts):
    I = counts.shape[0]
    idx_keep = []
    
    for i in xrange(0, I):
        a_N = counts[i, 0]
        b_N = counts[i, 1]
        d_N = a_N + b_N
        BAF_N = b_N/d_N
        
        if BAF_N >= BAF_N_MIN and BAF_N <= BAF_N_MAX:
            idx_keep.append(i)
            
    counts = counts[idx_keep]
    
    return counts

def tumor_LOH_test(counts):
    I = counts.shape[0]
    
    LOH_num = 0.0
    
    for i in xrange(0, I):
        a_T = counts[i, 2]
        b_T = counts[i, 3]
        d_T = a_T + b_T
        BAF_T = b_T/d_T
        
        if BAF_T < BAF_T_MIN or BAF_T > BAF_T_MAX:
            LOH_num = LOH_num + 1
    
    LOH_frec = LOH_num/I
    
    if LOH_frec < LOH_FREC_THRED:
        LOH_flag = False
    else:
        LOH_flag = True
        
    return (LOH_flag, LOH_frec)

def get_mu_H(mu_N, mu_T, phi):
    
    return (1 - phi)*mu_N + phi*mu_T
    
def get_phi(mu_N, mu_T, mu_H):
    
    return (mu_H - mu_N)/(mu_T - mu_N)
    
def get_phi_c(mu_N, mu_T, c_N, c_T, mu_H):
    
    return (mu_N - mu_H)*c_N/((mu_H - mu_T)*c_T + (mu_N - mu_H)*c_N)

def log_likelihoods(counts, phi, rho):
    I = counts.shape[0]
    G = GENOTYPES_TUMOR_NUM
    
    b_T = counts[:, 3]
    d_T = counts[:, 2] + counts[:, 3]
    
    mu_H = np.ones(G)
    
    for g in range(0, G):
        mu_H[g] = get_mu_H(MU_N[0], MU_T[g], phi)
        
    ll = np.log(rho) + log_binomial_likelihood(b_T, d_T, mu_H)
    
    return ll

def main():
    infile = sys.argv[1]
    iters = int(sys.argv[2])
    
################################################################################
#   IO and filtering
################################################################################
    
    reader = JointCountsReader(infile)
    
    chrom_set = reader.get_table_list()
    
    counts = []
    I = []
    segment_list = []
    
    for chrom in CHROM_LIST:
        if chrom not in chrom_set:
            continue
        
        counts_j = reader.get_counts(chrom)*1.0
        
        counts_j = normal_heterozygous_filter(counts_j)
        
        if counts_j.shape[0] < SITES_NUM_MIN:
            continue
        
        LOH_flag, LOH_frec = tumor_LOH_test(counts_j)
        
        if LOH_flag == True:
            counts.append(counts_j)
            I.append(counts_j.shape[0])
            segment_list.append(chrom)
            print '%s\t%s\t%s' % (chrom, LOH_frec, counts_j.shape[0]) 

################################################################################
#   Initialization
################################################################################

    J = len(counts)
    G = GENOTYPES_TUMOR_NUM
    phi = PHI_INIT
    rho = np.ones((J, G))*1.0/G
    I = np.array(I)

################################################################################
#   EM
################################################################################
    for t in range(0, iters):        
        if t <= BURN_IN:
            eta = ETA
        else:
            eta = 1
            
        #E-step
        u = np.zeros((J, G))
        v = np.zeros((J, G))
        w = np.zeros((J, G))
        
        for j in range(0, J):
            b_T_j = counts[j][:, 3]
            d_T_j = counts[j][:, 2] + counts[j][:, 3]
            
            log_likelihoods_j = log_likelihoods(counts[j], phi, rho[j])
            #xi_j = log_space_normalise_rows_annealing(log_likelihoods_j, eta)
            xi_j = log_space_normalise_rows(log_likelihoods_j)
            
            u_j = np.add.reduce(xi_j, axis = 0)
            v_j = xi_j*np.dot(b_T_j.reshape(I[j], 1), np.ones((1, G)))
            v_j = np.add.reduce(v_j, axis = 0)
            w_j = xi_j*np.dot(d_T_j.reshape(I[j], 1), np.ones((1, G)))
            w_j = np.add.reduce(w_j, axis = 0)
            
            u[j] = u_j
            v[j] = v_j
            w[j] = w_j
            
        #M-step
        phi = 0
        weights_sum = 0
        mu_H = np.zeros((J, G))
        phi_est = np.ones((J, G))*-1
        
        for j in range(0, J):
            rho[j] = (u[j] + omega - 1)/(I[j] + (omega - 1).sum())
            mu_H_j = v[j]/w[j]
            mu_H[j] = mu_H_j
            
            for g in range(0, G):
                if MU_N[0] == MU_T[g]:
                    continue
                
                if mu_H_j[g] < ERR:
                    mu_H_j[g] = ERR
                elif mu_H_j[g] > 1 - ERR:
                    mu_H_j[g] = 1 - ERR
                    
                #phi_j_g = get_phi(MU_N[0], MU_T[g], mu_H_j[g])
                phi_j_g = get_phi_c(MU_N[0], MU_T[g], C_N[0], C_T[g], mu_H_j[g])
                
                phi_est[j, g] = phi_j_g
                
                if phi_j_g < 0 or phi_j_g > 1:
                    continue
                    
                phi = phi + rho[j, g]*I[j]*phi_j_g
                weights_sum = weights_sum + rho[j, g]*I[j]
                
        phi = phi/weights_sum
        
        #log-likelihoods
        
        ll = 0
        for j in range(0, J):            
            log_likelihoods_j = log_likelihoods(counts[j], phi, rho[j])
            ll_j = np.logaddexp.reduce(log_likelihoods_j, axis = 1)
            ll = ll + np.add.reduce(ll_j, axis = 0)
        
        
        print '##########iter %s###########' % (t+1)
        print 'phi: %s\tll: %s' % (phi, ll)
    
    print 'rho:'
    print rho
    print 'mu_H:'
    print mu_H
    print 'phi_est'
    print phi_est
    print 'u'
    print u



    reader.close()

if __name__ == '__main__':
    main()