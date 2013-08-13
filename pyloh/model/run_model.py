'''
Created on 2012-08-13

@author: Yi Li
'''
from pyloh.model.poisson_model import PoissonProbabilisticModel

def run_model(args):
    model = PoissonProbabilisticModel()
    model.read_priors(args.priors_file_name)
    model.read_data(args.data_file_basename)
    model.preprocess_data()
    model.run(args.max_iters, args.stop_value)
    
