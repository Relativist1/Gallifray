# Saurabh
# mcmc embedding 1st test 4:36am 21st Jan, 2021
#---------------------------------------------------------------------------------#

import numpy as np
import matplotlib.pyplot as plt
import corner
import sys
import gallifray as gr
from gallifray.models import xsring
import ehtim as eh

#----------Variables that can be chaged----------#
I0_true = 2
Rp_true = 37
Rn_true = 22
ecn_true = 0.6
f_true = 0.3
phi_true = np.pi/2
fov = (Rp_true*2) + 20

n_walkers = 1000
n_samples =  5000

model_type = 'xsring'
likelihood = 'default'
#---------------------------------------------------------------------------------------
path = '/Users/geodesix/Desktop/Non-Kerr/2019-D01-01-master/uvfits/SR1_M87_2017_101_hi_hops_netcal_StokesI.uvfits'
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

#---------------------------------------------------------------------------------------

def ln_prob(param, uv, obs_amp, obs_sigma, model_type, model_fov, likelihood, prior_type):
    f = gr.mcmc_utils.ln_probability(param=param,
                                    uv=UV,
                                    obs_amp=vis_amp,
                                    obs_sigma=vis_err,
                                    model_type=model_type,
                                    model_fov=fov,
                                    likelihood=likelihood,
                                    prior_type='uniform',
                                    )
    return f

#---------------------------------------------------------------------------------------

prepare_mcmc = gr.mcmc(initial, likelihood=ln_prob, model_type=model_type, uv=UV,
                      obs=vis_amp, obs_sigma=vis_err, model_fov=fov,
                      n_walkers=n_walkers, n_samples=n_samples,
                      prior_type='uniform', likelihood_type='gaussian')
          
sampler = prepare_mcmc.run_sampler()
samples = sampler.get_chain()

ndim = len(initial)
fig, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True)

labels = [r"$I_0$", r"$R_p (\mu as)$", "$R_n (\mu as)$", r"${\epsilon}$", r"f", r"${\phi}$"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])

axes[-1].set_xlabel("step number");
fig.savefig('test1.png')
final = sampler.get_chain(flat=True)

#fig = corner.corner(final, labels=labels, truths=initial,show_titles=True)
#plt.show()
fig.savefig('cr_test1.png')

gr.Tardis(final, labels=labels, truths=initial, savefig='cr_grtest1.png', diag_shade_color='Red')
plt.show()
