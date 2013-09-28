'''
Created on 2010-08-30 for log_binomial_likelihood
Created on 2010-11-22 for log_space_normalise_rows_annealing

@author: Andrew Roth

JointSNVMix-0.6.2
joint_snv_mix.classification.utils.log_pdf.log_binomial_likelihood
joint_snv_mix.classification.utils.normalise.log_space_normalise_rows

================================================================================

Modified on 2013-07-31

@author: Yi Li
'''
import numpy as np
from scipy.special import gammaln

from pyloh import constants

#JointSNVMix
def log_space_normalise_rows_annealing(log_X, eta):
    nrows = log_X.shape[0]
    shape = ( nrows, 1 )
    log_X = log_X*eta
    
    log_norm_const = np.logaddexp.reduce( log_X, axis=1 )
    log_norm_const = log_norm_const.reshape( shape )

    log_X = log_X - log_norm_const
    
    X = np.exp( log_X )
    
    dt = X.dtype
    eps = np.finfo( dt ).eps
    
    X[X <= eps] = 0.
    
    return X

#JointSNVMix
def log_binomial_likelihood(k, n, mu):
    column_shape = (k.size, 1)
    k = k.reshape(column_shape)
    n = n.reshape(column_shape)
    
    row_shape = (1, mu.size)
    mu = mu.reshape(row_shape)
    
    return k * np.log(mu) + (n - k) * np.log(1 - mu)

def log_poisson_likelihood(k, Lambda):
    row_shape = (1, Lambda.size)
    Lambda = Lambda.reshape(row_shape)
    
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)

def get_genotypes_tumor(allele_number_max):
    genotypes_tumor = []
    
    for B_num in range(0, allele_number_max + 1):
        for A_num in range(0, allele_number_max + 1):
            if A_num == 0 and B_num == 0:
                g_T = 'NULL'
            else:
                g_T = 'A'*A_num + 'B'*B_num
                
            genotypes_tumor.append(g_T)
    
    return genotypes_tumor
    
def get_alleletypes_tumor(allele_number_max):
    alleletypes_tumor = []

    for B_num in range(0, allele_number_max + 1):
        for A_num in range(0, allele_number_max + 1):
            if A_num == 0 and B_num == 0:
                allele_T = 'NULL'
            elif A_num == B_num:
                allele_T = 'A'*A_num + 'B'*B_num
            elif A_num < B_num:
                allele_T = 'A'*B_num + 'B'*A_num + '/' + 'A'*A_num + 'B'*B_num
            else:
                allele_T = 'A'*A_num + 'B'*B_num + '/' + 'A'*B_num + 'B'*A_num
                
            alleletypes_tumor.append(allele_T)
    
    return alleletypes_tumor

def get_copynumber_tumor(allele_number_max):
    copynumber_tumor = []
    
    for B_num in range(0, allele_number_max + 1):
        for A_num in range(0, allele_number_max + 1):
            c_T = A_num + B_num
                
            copynumber_tumor.append(c_T)
    
    return copynumber_tumor

def get_MU_T(allele_number_max):
    empiri_BAF = constants.EMPIRI_BAF
    empiri_AAF = constants.EMPIRI_AAF
    err = constants.ERR
    
    MU_T = []
    
    for B_num in range(0, allele_number_max + 1):
        for A_num in range(0, allele_number_max + 1):
            if A_num == 0 and B_num == 0:
                mu_T = empiri_BAF/(empiri_BAF + empiri_AAF)
            elif A_num == 0 and B_num != 0:
                mu_T = 1 - err
            elif A_num != 0 and B_num == 0:
                mu_T = err
            else:
                mu_T = empiri_BAF*B_num/(empiri_BAF*B_num + empiri_AAF*A_num)
                
            MU_T.append(mu_T)
    
    return MU_T

def get_genotypes_tumor_num(allele_number_max):
    
    return (allele_number_max + 1)*(allele_number_max + 1)

def get_omega(allele_number_max):
    G = get_genotypes_tumor_num(allele_number_max)
    
    omega = [10 for i in range(0, G)]
    
    return omega

def get_x_E(x_N, x_T, phi):
    
    return (1 - phi)*x_N + phi*x_T
    
def get_c_E(c_N, c_T, phi):
    
    return (1 - phi)*c_N + phi*c_T

def get_mu_E(mu_N, mu_T, c_N, c_T, phi):
    
    return ((1 - phi)*c_N*mu_N + phi*c_T*mu_T)/((1 - phi)*c_N + phi*c_T)
    
def get_phi(x_N, x_T, x_E):
    
    return (x_E - x_N)/(x_T - x_N)
    
def get_phi_CNV(c_N, c_T, D_N, D_T, c_S, Lambda_S):
    numerator = c_N*D_T/(Lambda_S*D_N) - c_N
    denominator = c_T - c_N - c_S*D_T/(Lambda_S*D_N) + c_N*D_T/(Lambda_S*D_N)
    
    return numerator/denominator
    
def get_phi_LOH(mu_N, mu_T, mu_E, c_N, c_T):
    
    return (mu_N - mu_E)*c_N/((mu_E - mu_T)*c_T + (mu_N - mu_E)*c_N)