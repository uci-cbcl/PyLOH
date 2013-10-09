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

#JointSNVMix
def log_dirichlet_pdf(x, omega):
    if omega.ndim == 1: 
        log_normalisation_constant = gammaln(omega.sum()) - gammaln(omega).sum()
    
        log_likelihood = np.sum((omega - 1) * np.log(x))
    else:
        log_normalisation_constant = gammaln(omega.sum(axis=1)) - gammaln(omega).sum(axis=1)
    
        log_likelihood = np.sum((omega - 1) * np.log(x), axis=1)

    log_p = log_normalisation_constant + log_likelihood

    return log_p

def log_poisson_likelihood(k, Lambda):
    row_shape = (1, Lambda.size)
    Lambda = Lambda.reshape(row_shape)
    
    return k * np.log(Lambda) - Lambda - gammaln(k + 1)




#def get_copynumber_tumor(copynumber_max):
#    
#    return range(0, copynumber_max + 1)
#
#def get_copynumber_tumor_num(copynumber_max):
#    
#    return copynumber_max + 1

#def get_genotypes_tumor(copynumber_max):
#    genotypes_tumor = []
#    
#    C_T = get_copynumber_tumor(copynumber_max)
#    
#    for c in C_T:
#        for B_num in range(0, c+1):
#            A_num = c - B_num
#            
#            if A_num == 0 and B_num == 0:
#                g_T = 'NULL'
#            else:
#                g_T = 'A'*A_num + 'B'*B_num
#                
#            genotypes_tumor.append(g_T)
#    
#    return genotypes_tumor
#
#def get_genotypes_tumor_num(copynumber_max):
#    
#    return len(get_genotypes_tumor(copynumber_max))

#def get_copynumber_tumor_compat(copynumber_max):
#    copynumber_tumor_compat = []
#    
#    C_T = get_copynumber_tumor(copynumber_max)
#    
#    for c in C_T:
#        for B_num in range(0, c+1):
#            A_num = c - B_num
#            
#            c_gT = A_num + B_num
#                
#            copynumber_tumor_compat.append(c_gT)
#    
#    return copynumber_tumor_compat

#def get_MU_T(copynumber_max):
#    empiri_BAF = constants.EMPIRI_BAF
#    empiri_AAF = constants.EMPIRI_AAF
#    err = constants.ERR
#    
#    MU_T = []
#    C_T = get_copynumber_tumor(copynumber_max)
#    
#    for c in C_T:
#        for B_num in range(0, c+1):
#            A_num = c - B_num
#            
#            if A_num == 0 and B_num == 0:
#                mu_T = empiri_BAF/(empiri_BAF + empiri_AAF)
#            elif A_num == 0 and B_num != 0:
#                mu_T = 1 - err
#            elif A_num != 0 and B_num == 0:
#                mu_T = err
#            else:
#                mu_T = empiri_BAF*B_num/(empiri_BAF*B_num + empiri_AAF*A_num)
#                
#            MU_T.append(mu_T)
#    
#    return MU_T

#def get_P_CG(copynumber_max):
#    sigma = constants.SIGMA
#    
#    C = get_copynumber_tumor_num(copynumber_max)
#    G = get_genotypes_tumor_num(copynumber_max)
#    c_T = np.array(get_copynumber_tumor(copynumber_max))
#    g_T = get_genotypes_tumor(copynumber_max)
#    
#    P_CG = np.ones((C, G))*sigma
#    
#    for c in range(0, C):
#        if c_T[c] == 0:
#            P_CG[c, 0] = 1 - sigma*(G - 1)
#            continue
#    
#        if c_T[c]%2 == 0:
#            compatible_num = 1
#            for g in range(0, G):
#                A_num = g_T[g].count('A')
#                B_num = g_T[g].count('B')
#                if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL' and A_num == B_num:
#                    P_CG[c, g] = (1 - sigma*(G - compatible_num))/compatible_num
#            
#            continue
#        
#        compatible_num = 0
#        for g in range(0, G):
#            if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL':
#                compatible_num += 1
#                        
#        for g in range(0, G):
#            if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL':
#                P_CG[c, g] = (1 - sigma*(G - compatible_num))/compatible_num
#    
#    return P_CG

#def get_omega(copynumber_max):
#    C = get_copynumber_tumor_num(copynumber_max)
#    
#    omega = [10 for i in range(0, C)]
#    
#    return omega

def get_genotypes_tumor(allelenumber_max):
    genotypes_tumor = []
    
    for B_num in range(0, allelenumber_max + 1):
        for A_num in range(0, allelenumber_max + 1):
            if A_num == 0 and B_num == 0:
                g_T = 'NULL'
            else:
                g_T = 'A'*A_num + 'B'*B_num
                
            genotypes_tumor.append(g_T)
    
    return genotypes_tumor

def get_genotypes_tumor_num(allelenumber_max):
    
    return (allelenumber_max + 1)*(allelenumber_max + 1)

def get_copynumber_tumor_compat(allelenumber_max):
    copynumber_tumor_compat = []
    
    for B_num in range(0, allelenumber_max + 1):
        for A_num in range(0, allelenumber_max + 1):
            c_gT = A_num + B_num
                
            copynumber_tumor_compat.append(c_gT)
    
    return copynumber_tumor_compat

def get_copynumber_tumor(allelenumber_max):
    
    return range(0, allelenumber_max*2 + 1)
    
def get_copynumber_tumor_num(allelenumber_max):
    
    return allelenumber_max*2 + 1

def get_MU_T(allelenumber_max):
    empiri_BAF = constants.EMPIRI_BAF
    empiri_AAF = constants.EMPIRI_AAF
    err = constants.ERR
    
    MU_T = []
    
    for B_num in range(0, allelenumber_max + 1):
        for A_num in range(0, allelenumber_max + 1):
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

def get_P_CG(allelenumber_max):
    sigma = constants.SIGMA
    
    C = get_copynumber_tumor_num(allelenumber_max)
    G = get_genotypes_tumor_num(allelenumber_max)
    c_T = np.array(get_copynumber_tumor(allelenumber_max))
    g_T = get_genotypes_tumor(allelenumber_max)
    
    P_CG = np.ones((C, G))*sigma
    
    for c in range(0, C):
        if c_T[c] == 0:
            P_CG[c, 0] = 1 - sigma*(G - 1)
            continue
        
        if c_T[c]%2 == 0:
            compatible_num = 1
            for g in range(0, G):
                A_num = g_T[g].count('A')
                B_num = g_T[g].count('B')
                if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL' and A_num == B_num:
                    P_CG[c, g] = (1 - sigma*(G - compatible_num))/compatible_num
            
            continue
        
        compatible_num = 0
        for g in range(0, G):
            if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL':
                compatible_num += 1
                
        for g in range(0, G):
            if len(g_T[g]) == c_T[c] and g_T[g] != 'NULL':
                P_CG[c, g] = (1 - sigma*(G - compatible_num))/compatible_num
    
    return P_CG

def get_omega(allelenumber_max):
    C = get_copynumber_tumor_num(allelenumber_max)
    
    omega = [10 for i in range(0, C)]
    
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
    
def get_phi_range_2nd(phi):
    
    return np.arange(phi - 0.01, phi + 0.01, 0.001).tolist()