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
