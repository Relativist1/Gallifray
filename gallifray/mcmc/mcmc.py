#    MCMC Class
#
#    Copyright (C) 2020 Saurabh
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#edit 9th April, 2022

from __future__ import division
from __future__ import print_function

from builtins import object
import sys, os
import numpy as np
from scipy.optimize import minimize
import emcee

import gallifray as gr

def blockp():
    sys.stdout = open(os.devnull,'w')

class mcmc(object):
    """MCMC sampling class
    
    Attributes:
        param (list) : list of parameters
        
    """
    def __init__(self, param, model_type, obs_data, model_fov, theta_G, n_walkers=10, n_samples=1000, blob_width = 1e-2, likelihood_type='gaussian', prior_type='uniform',exec_c='./grrt', interp=None):
        """
        
        Args:
            param (list): list of parameters
            
        Return:
            
        """

        u =  obs_data.unpack(['u'],conj=True)['u']
        v =  obs_data.unpack(['v'],conj=True)['v']
        obs = obs_data.unpack(['amp'],conj=True)['amp']
        obs_sigma = obs_data.unpack(['sigma'],conj=True)['sigma']
        uv = [u, v]

        self.obs_data = obs_data
        self.param = param
        self.model_type = model_type
        self.uv = uv
        self.obs = obs
        self.obs_sigma = obs_sigma
        self.likelihood_type = likelihood_type
        self.n_walkers = n_walkers
        self.n_samples = n_samples
        self.prior_type = prior_type
        self.model_fov = model_fov
        self.exec_c = exec_c
        self.interp = interp
        self.blob_width = blob_width
        self.theta_G = theta_G


    def run_sampler(self, n_threads="4", minimize_lk = False,mult_proc=False, progress=True, moves=None, args=None, kwargs=None, backend=None, vectorize=False, blobs_dtype=None, a=None, postargs=None, threads=None, **kkwargs):
        """An ensemble MCMC sample and samples iterator
        
        Args:
            mult_proc (bool): multiprocessing via paralleization, if True, then the sampler runs parallely.
        Return:
            
        """
        #edit 4th June, 2022
        init1 = []
        for i in self.param.keys():
            init1.append(self.param[i][0])

        initial = np.asarray(init1)
        ndim = len(initial)
        #edit 9th April, 2022

        if minimize_lk:

            def lk_h(p0, pr_param, exec_c):
                lk = gr.likelihood(p0, self.obs_data, self.model_type, self.model_fov, self.theta_G)
                return -lk.ln_gaussian_physical(pr_param, exec_c)

            sol1 = minimize(lk_h, initial, args=(self.param, self.exec_c))

            positions = sol1.x + self.blob_width*np.random.randn(self.n_walkers, ndim)
        else:

            positions = initial + self.blob_width*np.random.randn(self.n_walkers, ndim)

        # def ln_prob(p0, obs_data, model_type, model_fov,
        #     likelihood, prior_type,interp=None):
        #     f = ln_probability(p0=p0, pr_params=self.param, obs_data=obs_data,
        #                                     model=model_type,model_fov=self.model_fov,likelihood=likelihood,
        #                                     prior_type=prior_type,exec_c=self.exec_c, theta_G=self.theta_G, interp=interp)
        #     return f

        if mult_proc:

            from multiprocessing import Pool
            os.environ["OMP_NUM_THREADS"] = n_threads

            with Pool() as pool:
                samp = emcee.EnsembleSampler(self.n_walkers, ndim, ln_probability,pool=pool,
                                            args=(self.param, self.obs_data,
                                                self.model_type, self.model_fov, self.likelihood_type,
                                                self.prior_type, self.theta_G, self.exec_c),
                                            moves=moves, kwargs=kwargs, backend=backend,
                                            vectorize=vectorize, blobs_dtype=blobs_dtype,
                                            a=a, postargs=postargs, threads=threads)
                
                samp.run_mcmc(positions, self.n_samples, progress=progress, **kkwargs)
        else:
            samp = emcee.EnsembleSampler(self.n_walkers, ndim, ln_probability,
                                            args=(self.param, self.obs_data,
                                                self.model_type, self.model_fov, self.likelihood_type,
                                                self.prior_type, self.theta_G, self.exec_c),
                                            moves=moves, kwargs=kwargs, backend=backend,
                                            vectorize=vectorize, blobs_dtype=blobs_dtype,
                                            a=a, postargs=postargs, threads=threads)
        
            samp.run_mcmc(positions, self.n_samples, progress=progress, **kkwargs)
        
        
        
        return samp

def priori_transform(param, model, p_type='uniform'):
    pr = gr.priors(param)
    prr = pr.lnprior(model, p_type)
    return prr



def ln_probability(p0, pr_params, obs_data, model, model_fov, likelihood, prior_type, theta_G, exec_c='./grrt', interp=None):
    """Defines the log probability

    Return:
       log likelihood
    """
    def lkh(p0, pr_params, obs_data, model, model_fov, likelihood, interp=None):
        lk = gr.likelihood(p0, obs_data, model, model_fov, theta_G)
        if model[0]=='physical':
            lkhood = lk.ln_gaussian_physical(pr_param=pr_params, exec_c=exec_c)
        else:
            lkhood = lk.ln_gaussian(interp=interp)

         # likelihood=='rice' or likelihood=='ln_vis_amp':
         #    lkhood = lk.ln_vis_amp(interp=interp)
        return lkhood

    
    lnp = priori_transform(p0, model, prior_type)
    if not np.isfinite(lnp):
       return -np.inf
    return lnp + lkh(p0, pr_params, obs_data, model, model_fov, likelihood, interp)
