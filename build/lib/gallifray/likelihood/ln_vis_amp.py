import numpy as np
from gallifray.models import *
from scipy.special import i0

def vis_amp(param,obs_bl, obs_amp, obs_sigma, model_type):
    """GALLIFRAY :: Rice distribution likelihood for visibility amplitudes

    Args:
        param (list) : List of parameters
        obs_bl (array) : observation baselines (in uas)
        uv (list) : List of observation u,v points
        obs_sigma (array) : observation errors
        model_type (str) : Type of image model

    Return:
        Calculated Likelihood Distribution
    """

    if model_type=='sym_gauss':
        I0, sigma, fov = param
        dim = len(obs_amp)
        imarr = gauss(I0, sigma, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, interp=None)
        model_amp = model_vis['amp']
        
    if model_type=='asym_gauss':
        I0, A, sigma, phi, fov = param
        dim = len(obs_amp)
        imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, interp=None)
        model_amp = model_vis['amp']
        
    if model_type=='disk':
        I0, R, fov = param
        dim = len(obs_amp)
        imarr = disk(I0, R, fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, interp=None)
        model_amp = model_vis['amp']
        
    if model_type=='xsring':
        I0, R_p, R_n, ecn, f, phi, fov = param
        dim = len(obs_amp)
        imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, interp=None)
        model_amp = model_vis['amp']
        
    if model_type=='xsringauss':
        I0, R_p, R_n, ecn, f, gax, aq, gq, phi, fov = param
        dim = len(obs_amp)
        imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
        model_vis = imarr.vis_data(fov=fov, uv=uv, interp=None)
        model_amp = model_vis['amp']


    model = model_amp
    obs = obs_mp
    sigma = obs_sigma
    
    z = model * obs/sigma**2)
    lp = np.sum(np.log(model) - 2*np.log(sigma) -  (model + obs)**2/(2*sigma**2) + np.log(i0(z))
    return lp
