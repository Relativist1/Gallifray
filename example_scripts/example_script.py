#!/usr/bin/env python

"""example_script.py: Esitmation of best fit parameters using MCMC sampling of
a geometric model over supplied some observed dataset.
   
 ------------------------------- How to use --------------------------------
 run this as a usual python script in the shell.
 
 $ python example_script.py
    
 
 
 """

__author__      = "Barack Obama"
__copyright__   = "Copyright 2009, Planet Earth"

#---------------------------------------------------------------------------------#
import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
import gallifray as gr
from gallifray.models import xsring
import ehtim as eh


path = '../datasets/xsringauss_sgrA.uvfits' # path of the dataset
save_fig_dist = 'cr_grtest1.png'            # filename for saving distribution plot

#------------------------Inital guess for the parameters--------------------------#

I0_true = 2   # Refer to the paper for detailed explaination of model and parameters.
Rp_true = 37
Rn_true = 22
ecn_true = 0.6
f_true = 0.3
phi_true = np.pi/2
fov = (Rp_true*2) + 20

n_walkers = 500           # No of walkers
n_samples = 1000          # No of iterations

model_type = 'xsring'     # Define the model type
likelihood = 'ln_vis_amp' # Define the likelihood
prior_type = 'uniform'    # Define the type of prior

#--------------------------- Unpacking the dataset --------------------------------#

obs_m = eh.obsdata.load_uvfits(path)
obs_m.add_scans()
obs_m = obs_m.avg_coherent(0.,scan_avg=True)

u_eh =  obs_m.unpack(['u'],conj=True)['u']
v_eh =  obs_m.unpack(['v'],conj=True)['v']
vis_amp = obs_m.unpack(['amp'],conj=True)['amp']
vis_err = obs_m.unpack(['sigma'],conj=True)['sigma']
obs_bl = obs_m.unpack(['uvdist'],conj=True)['uvdist']
    
UV = [u_eh, v_eh]
initial = [I0_true, Rp_true, Rn_true, ecn_true, f_true, phi_true]

#----------------------------Define the likeliihod --------------------------------#
def ln_prob(param, uv, obs_amp, obs_sigma, model_type, model_fov,
            likelihood, prior_type):
    f = gr.mcmc_utils.ln_probability(param=param,
                                    uv=UV,
                                    obs_amp=vis_amp,
                                    obs_sigma=vis_err,
                                    model_type=model_type,
                                    model_fov=fov,
                                    likelihood=likelihood,
                                    prior_type=prior_type,
                                    )
    return f
    
#---------------------------Intialising the MCMC sampler---------------------------#

prepare_mcmc = gr.mcmc(initial, likelihood=ln_prob, model_type=model_type, uv=UV,
                      obs=vis_amp, obs_sigma=vis_err, model_fov=fov,
                      n_walkers=n_walkers, n_samples=n_samples,
                      prior_type=prior_type, likelihood_type=likelihood_type)
          
sampler = prepare_mcmc.run_sampler()
ndim = len(initial)
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)

labels = [r"$I_0$", r"$R_p (\mu as)$", "$R_n (\mu as)$", r"${\epsilon}$", r"f", r"${\phi}$"]

post = sampler.get_chain(flat=True)

#-----------Plot the Triangular Posterior-Distribution Plot (Takes a while)--------#

gr.Tardis(post, labels=labels, truths=initial,
         savefig=save_fig_dist,
         diag_shade_color='Red')
plt.show()
