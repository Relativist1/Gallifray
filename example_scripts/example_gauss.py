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
path = 'synthetic_data_asym_gauss.uvfits' # path of the dataset

n_walkers = 20           # No of walkers
n_samples = 500          # No of iterations

#------------------------Inital guess for the parameters--------------------------#
I0_true = 2 
sigma_true = 20
A_true = 0.5
phi_true = np.pi/2
fov = sigma_true*2 + 40


param = {'I0': [I0_true, 0, 10, 0.2],
         'sigma': [sigma_true, 0, 100, 0.3],
         'A' : [A_true, 0, 1, 0.2],
         'phi': [phi_true, -np.pi, np.pi, 0.2]
}

model_type =['geom','asym_gauss', param]   # Define the model type
likelihood = 'gaussian'         # Define the likelihood
prior_type = 'uniform'          # Define the type of prior

#--------------------------- Unpacking the dataset --------------------------------#

obs_m = eh.obsdata.load_uvfits(path)
obs_m.add_scans()
obs_m = obs_m.avg_coherent(0.0,scan_avg=True)

#---------------------------Intialising the MCMC sampler---------------------------#
filename = "chain_test_asym_gauss.h5"
backend = emcee.backends.HDFBackend(filename)
prepare_mcmc = gr.mcmc(param,fov_m=40, model_type=model_type,obs_data=obs_m, model_fov=fov,
                      n_walkers=n_walkers, n_samples=n_samples,
                      prior_type=prior_type, likelihood_type=likelihood, use_priori_default=True)
          
sampler = prepare_mcmc.run_sampler(init_position='guess',backend=backend)

# ndim = len(initial)
# # labels = [r"$I_0$", r"$R_p (\mu as)$", "$R_n (\mu as)$", r"${\epsilon}$", r"f", r"${\phi}$"]
# labels = [r"$I_0$", r"$\A (\mu as)$"]
# post = sampler.get_chain(flat=True)

# #-----------Plot the Triangular Posterior-Distribution Plot (Takes a while)--------#

# # gr.Tardis(post, labels=labels, truths=initial,
# #          savefig=save_fig_dist,
# #          diag_shade_color='Red')
# fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)

# samples = sampler.get_chain()

# for i in range(ndim):
#     ax = axes[i]
#     ax.plot(samples[:, :, i], "k", alpha=0.3)
#     ax.set_xlim(0, len(samples))
#     ax.set_ylabel(labels[i])
#     ax.yaxis.set_label_coords(-0.1, 0.5)

# axes[-1].set_xlabel("step number");

# plt.show()

# # tau = sampler.get_autocorr_time()
# # print(tau)

# final = sampler.get_chain(flat=True)
# #final = sampler.get_chain(discard=100, thin=20, flat=True)

# fig = corner.corner(final, labels=labels, truths=initial,show_titles=True);
# plt.show()
