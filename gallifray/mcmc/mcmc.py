from __future__ import division
from __future__ import print_function

from builtins import object
import sys, os

import numpy as np
from scipy.optimize import minimize
import emcee

import gallifray as gr
from gallifray.mcmc.mcmc_utils import *

class mcmc(object):
    """MCMC sampling class
    
    Attributes:
        param (list) : list of parameters
        
    """
    def __init__(self, param, likelihood, model_type, uv, obs, obs_sigma, model_fov, n_walkers=500, n_samples=1000, likelihood_type='gaussian', prior_type='uniform'):
        self.param = param
        """
        
        Args:
            param (list): list of parameters
            
        Return:
            
        """
        self.param = param
        self.model_type = model_type
        self.uv = uv
        self.obs = obs
        self.obs_sigma = obs_sigma
        self.likelihood = likelihood
        self.likelihood_type = likelihood_type
        self.n_walkers = n_walkers
        self.n_samples = n_samples
        self.prior_type = prior_type
        self.model_fov = model_fov


    def run_sampler(self, n_threads="1",mult_proc=False, progress=True, moves=None, args=None, kwargs=None, backend=None, vectorize=False, blobs_dtype=None, a=None, postargs=None, threads=None, **kkwargs):
        """An ensemble MCMC sample and samples iterator
        
        Args:
            mult_proc (bool): multiprocessing via paralleization, if True, then the sampler runs parallely.
        Return:
            
        """

        initial = np.asarray(self.param)
        ndim = len(initial)
        positions = initial + 0.5*np.random.randn(self.n_walkers, ndim)
        
        
        samp = emcee.EnsembleSampler(self.n_walkers, ndim, self.likelihood,
                                    args=(self.uv, self.obs, self.obs_sigma,
                                    self.model_type, self.model_fov,
                                    self.likelihood_type, self.prior_type),
                                    moves=moves, kwargs=kwargs, backend=backend,
                                    vectorize=vectorize, blobs_dtype=blobs_dtype,
                                    a=a, postargs=postargs, threads=threads)
        
        final = samp.run_mcmc(positions, self.n_samples, progress=progress, **kkwargs)
        
        
        if mult_proc==True:
            from multiprocessing import Pool
            os.environ["OMP_NUM_THREADS"] = n_threads
            
            if __name__ == '__main__':
                with Pool() as pool:
                    samp = emcee.EnsembleSampler(self.n_walkers, ndim,
                                                self.likelihood,pool=pool,
                                                args=(self.uv, self.obs, self.obs_sigma,
                                                self.model_fov, self.likelihood_type,
                                                self.prior_type),
                                                moves=moves, kwargs=kwargs, backend=backend,
                                                vectorize=vectorize, blobs_dtype=blobs_dtype,
                                                a=a, postargs=postargs, threads=threads)
                    
                    samp.run_mcmc(positions, self.n_samples, progress=progress, **kkwargs)
                
        return samp
