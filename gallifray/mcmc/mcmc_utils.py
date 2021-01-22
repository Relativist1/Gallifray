#    MCMC utility functions
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
import numpy as np
import gallifray as gr

def priori_transform(param, model_type, type='uniform'):
    pr = gr.priors(param)
    prr = pr.lnprior(model_type, type=type)
    return prr




def ln_probability(param, uv, obs_amp, obs_sigma, model_type, model_fov, likelihood, prior_type):
    """Defines the log probability

    Return:
       log likelihood
    """
    def lkh(param, uv, obs_amp, obs_sigma, model_type, model_fov, likelihood='default'):
        lk = gr.likelihood(param, uv, obs_amp, obs_sigma, model_type, model_fov)
        lkhood = lk.ln_gaussian()
        return lkhood
    
    lnp = priori_transform(param, model_type, prior_type)
    if not np.isfinite(lnp):
       return -np.inf
    return lnp + lkh(param, uv, obs_amp, obs_sigma, model_type, model_fov, likelihood)
