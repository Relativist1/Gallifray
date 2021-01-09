import numpy as np
from gallifray.models import *
from gallifray.utilities.random import model_datapoints

def gaussian(param,obs_bl, uv, obs_amp, obs_sigma, model_type, **kwargs):
    """GALLIFRAY :: Gaussian distribution likelihood for general purposes
    
    Args:
        param (list) : List of parameters
        obs_bl (array) : observation baselines (in uas)
        uv (list) : List of observation u,v points
        obs_sigma (array) : observation errors
        model_type (str) : Type of image model
        **kwargs : other arguments for the visibility model object
    Return:
        Calculated Likelihood Distribution
    """
    
    if model_type=='sym_gauss':
        I0, sigma, fov = param
        dim = len(obs_amp)
        imarr = gauss(I0, sigma, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, **kwargs)
        model_amp = model_vis['amp']
        
    if model_type=='asym_gauss':
        I0, A, sigma, phi, fov = param
        dim = len(obs_amp)
        imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, **kwargse)
        model_amp = model_vis['amp']
        
    if model_type=='disk':
        I0, R, fov = param
        dim = len(obs_amp)
        imarr = disk(I0, R, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, **kwargs)
        model_amp = model_vis['amp']
        
    if model_type=='xsring':
        I0, R_p, R_n, ecn, f, phi, fov = param
        dim = len(obs_amp)
        imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, **kwargs)
        model_amp = model_vis['amp']
        
    if model_type=='xsringauss':
        I0, R_p, R_n, ecn, f, gax, aq, gq, phi, fov = param
        dim = len(obs_amp)
        imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, **kwargs)
        model_amp = model_vis['amp']

    model_new_amp = model_datapoints(index=Index, obs_uvdist=obs_bl, model_amp=model_amp)
    model = model_new_amp
    obs = obs_new_amp
    sigma = obs_new_sigma
    lp = -0.5*sum((model - obs)**2/sigma + np.log(sigma) + np.log(2*np.pi))
    return lp
