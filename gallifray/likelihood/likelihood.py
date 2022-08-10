#    Likelihood Class
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

from __future__ import division
from __future__ import print_function

from builtins import object

import numpy as np
from scipy.special import i0
import os
import subprocess
from gallifray.models import *
import ehtim as eh
import sys

def blockp():
    sys.stdout = open(os.devnull,'w')

class likelihood(object):
    """MCMC sampling class
    
    Attributes:
        param (list) : list of parameters
        model (str) : model type
        
    """
    def __init__(self, param, obs_data, model, model_fov, theta_G):
        """
        
        Args:
            param (list) : List of parameters
            uv (list) : List of observation u,v points
            obs_sigma (array) : Observation errors
            model (str) : Type of image model
            
        Return:
            
        """
        u =  obs_data.unpack(['u'],conj=True)['u']
        v =  obs_data.unpack(['v'],conj=True)['v']
        obs_amp = obs_data.unpack(['amp'],conj=True)['amp']
        obs_sigma = obs_data.unpack(['sigma'],conj=True)['sigma']
        uv = [u, v]

        self.param = param
        self.uv = uv
        self.obs_amp = obs_amp
        self.obs_sigma = obs_sigma
        self.model = model
        self.model_fov = model_fov
        self.theta_G = theta_G
        self.obs_data = obs_data

    def ln_gaussian(self, interp, **kwargs):
        """GALLIFRAY :: Gaussian distribution likelihood for general purposes
        
        Args:
            **kwargs : other arguments for the visibility model object
        Return:
            Calculated Likelihood Distribution
        """
        
        fov = self.model_fov
        if self.model=='sym_gauss':
            I0, sigma = self.param
            dim = len(self.obs_amp)
            imarr = gauss(I0, sigma, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='asym_gauss':
            I0, A, sigma, phi = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='disk':
            I0, R = self.param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='crescent':
            I0, R, psi, tau, phi = self.param
            dim = len(self.obs_amp)
            imarr = crescent(I0, R, psi, tau, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='xsring':
            I0, R_p, R_n, ecn, f, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        model = model_amp
        obs = self.obs_amp
        sigma = self.obs_sigma
        
        lp = -0.5*sum((model - obs)**2/sigma + np.log(sigma) + np.log(2*np.pi))
        
        return lp
        
    def ln_vis_amp(self, interp, **kwargs):
        """GALLIFRAY :: Rice distribution likelihood for visibility amplitudes

        Args:
            **kwargs : other arguments for the visibility model object
        Return:
            Calculated Likelihood Distribution
        """
        fov = self.model_fov
        if self.model=='sym_gauss':
            I0, sigma = self.param
            dim = len(self.obs_amp)
            imarr = gauss(I0, sigma, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='asym_gauss':
            I0, A, sigma, phi = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargse)
            model_amp = model_vis['amp']
            
        if self.model=='disk':
            I0, R = self.param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='crescent':
            I0, R, psi, tau, phi = self.param
            dim = len(self.obs_amp)
            imarr = crescent(I0, R, psi, tau, phi, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='xsring':
            I0, R_p, R_n, ecn, f, phi  = self.param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, interp=interp, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']


        model = model_amp
        obs = self.obs_amp
        sigma = self.obs_sigma
        z = (model * obs)/sigma**2
        
        lp = sum(np.log(model) - 2*np.log(sigma) -  (model + obs)**2/(2*sigma**2) + np.log(i0(z)))/(dim + len(self.param))
        return lp

    def ln_gaussian_physical(self, pr_param, exec_c, fov_m):
        """GALLIFRAY :: Gaussian distribution likelihood for general purposes
        
        Args:
            exec_c (list): executable along with extra arguments
        Return:
            Calculated Likelihood Distribution
        """

        X, Y, Z = prep_grrt(self.param, exec_c, pr_param, fov_m) 

        model_amp, err = vis_observe(X, Y, Z, self.obs_data, self.theta_G)

        model = model_amp
        obs = self.obs_amp
        sigma = self.obs_sigma
        
        lp = -0.5*sum((model - obs)**2/sigma + np.log(sigma) + np.log(2*np.pi))

        print(lp)
        return lp

def prep_grrt(p0, exec_c, fov_m):

    get_args = exec_c

    #decompose and creates a list of all arguments
    for i in range(len(p0)):
        get_args.append("{}".format(p0[i]))
    
    out_rt = subprocess.Popen(get_args, stdout=subprocess.PIPE).communicate()
    # print("done\n")
    X = []
    Y = []
    Z = []   
    

    for i in out_rt[0].decode('utf-8').splitlines():
        X.append(float(i.split(",")[0]))
        Y.append(float(i.split(",")[1]))
        Z.append(float(i.split(",")[2]))
    
    x_im = np.asarray(X)
    y_im = np.asarray(Y)
    z_im = np.asarray(Z)

    dim = int(np.sqrt(len(X)))
    

    x1  = fov_m/2;
    x2  = fov_m/(dim+1.)
    x3 = -x1 + x2*(x_im)
    y3 = -x1 + x2*(y_im)
    
    x_im = x3.reshape((dim, dim))
    y_im = y3.reshape((dim, dim))
    z_im = z_im.reshape((dim, dim))

    return [x_im, y_im,  z_im]

def vis_observe(X, Y, Z, obs_ref, theta_G):
    ra = obs_ref.ra
    dec = obs_ref.dec
    rf = obs_ref.rf
    source = obs_ref.source 

    psize = abs(Y[0][0] - Y[0][1])*eh.RADPERUAS*theta_G

    im1 = eh.image.Image(np.flipud(Z),psize=psize,ra=ra, dec=dec, rf= rf, source=source)

    obs = im1.observe_same(obs_ref, add_th_noise=True,verbose=False)
    blockp()
    vis_amp1 = obs_ref.unpack(['amp'],conj=True)['amp']
    vis_err1 = obs_ref.unpack(['sigma'],conj=True)['sigma'] 
    
    return vis_amp1, vis_err1