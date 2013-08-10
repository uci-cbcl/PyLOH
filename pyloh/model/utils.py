'''
Created on 2013-07-31

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

def get_x_E(x_N, x_T, phi):
    
    return (1 - phi)*x_N + phi*x_T
    
def get_phi(x_N, x_T, x_E):
    
    return (x_E - x_N)/(x_T - x_N)
    
def get_phi_CNV(c_N, c_T, c_E):
    
    return (c_E - c_N)/(c_T - c_N)
    
def get_phi_LOH(mu_N, mu_T, mu_E, c_N, c_T):
    
    return (mu_N - mu_E)*c_N/((mu_E - mu_T)*c_T + (mu_N - mu_E)*c_N)