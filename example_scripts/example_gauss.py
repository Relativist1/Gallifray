#!/usr/bin/env python

"""example_script.py: Esitmation of best fit parameters using MCMC sampling of
a geometric model over supplied some observed dataset.
   
 ------------------------------- How to use --------------------------------
 run this as a usual python script in the shell.
 
 $ python example_script.py
    
 
 
 """

__author__      = "Saurabh"
__copyright__   = "Copyright 2020, Gallifray"

#---------------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
import gallifray as gr
from gallifray.models import xsring
import ehtim as eh
import corner
import emcee
path = 'test_asym_gauss.uvfits' # path of the dataset

n_walkers = 50           # No of walkers
n_samples = 500          # No of iterations

#------------------------Inital guess for the parameters--------------------------#
I0_true = 2 
sigma_true = 40
A_true = 0.5
phi_true = np.pi/2
fov = sigma_true*2 + 400


param = {'I0': [I0_true, 0, 10, 0.1],
         'sigma': [sigma_true, 0, 100, 0.3],
         'A' : [A_true, 0, 1, 0.2],
         'phi': [phi_true, 0, np.pi, 0.2]
}

model_type =['geom','asym_gauss', param]   # Define the model type
likelihood = 'gaussian'         # Define the likelihood
prior_type = 'uniform'          # Define the type of prior

#--------------------------- Unpacking the dataset --------------------------------#

obs_m = eh.obsdata.load_uvfits(path)
obs_m.add_scans()
obs_m = obs_m.avg_coherent(0.0,scan_avg=True)

#---------------------------Intialising the MCMC sampler---------------------------#
filename = "chain_test_asym_gauss_sim1.h5"
backend = emcee.backends.HDFBackend(filename)
prepare_mcmc = gr.mcmc(param,fov_m=40, model_type=model_type,obs_data=obs_m, model_fov=fov,
                      n_walkers=n_walkers, n_samples=n_samples,
                      prior_type=prior_type, likelihood_type=likelihood, use_priori_default=False)
          
sampler = prepare_mcmc.run_sampler(init_position='guess',backend=backend)