'''
Created on 2013-07-31

@author: Yi Li
'''
import numpy as np

from pyloh import constants

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