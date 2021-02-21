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

from gallifray.models import *

class likelihood(object):
    """MCMC sampling class
    
    Attributes:
        param (list) : list of parameters
        model_type (str) : model type
        
    """
    def __init__(self, param, uv, obs_amp, obs_sigma, model_type, model_fov):
        """
        
        Args:
            param (list) : List of parameters
            uv (list) : List of observation u,v points
            obs_sigma (array) : Observation errors
            model_type (str) : Type of image model
            
        Return:
            
        """
        self.param = param
        self.uv = uv
        self.obs_amp = obs_amp
        self.obs_sigma = obs_sigma
        self.model_type = model_type
        self.model_fov = model_fov
    def ln_gaussian(self, **kwargs):
        """GALLIFRAY :: Gaussian distribution likelihood for general purposes
        
        Args:
            **kwargs : other arguments for the visibility model object
        Return:
            Calculated Likelihood Distribution
        """
        
        fov = self.model_fov
        if self.model_type=='sym_gauss':
            I0, sigma = self.param
            dim = len(self.obs_amp)
            imarr = gauss(I0, sigma, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='asym_gauss':
            I0, A, sigma, phi = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='disk':
            I0, R = param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='xsring':
            I0, R_p, R_n, ecn, f, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']

        model = model_amp
        obs = self.obs_amp
        sigma = self.obs_sigma
        lp = -0.5*sum((model - obs)**2/sigma + np.log(sigma) + np.log(2*np.pi))
        return lp
        
    def ln_vis_amp(self, **kwargs):
        """GALLIFRAY :: Rice distribution likelihood for visibility amplitudes

        Args:
            **kwargs : other arguments for the visibility model object
        Return:
            Calculated Likelihood Distribution
        """
    
        if self.model_type=='sym_gauss':
            I0, sigma, fov = self.param
            dim = len(self.obs_amp)
            imarr = gauss(I0, sigma, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='asym_gauss':
            I0, A, sigma, phi, fov = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargse)
            model_amp = model_vis['amp']
            
        if self.model_type=='disk':
            I0, R, fov = param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='xsring':
            I0, R_p, R_n, ecn, f, phi, fov = param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
            
        if self.model_type=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi, fov = param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=dim)
            model_vis = imarr.vis_data(fov=fov, uv=self.uv, **kwargs)
            model_amp = model_vis['amp']
        
        model = model_amp
        obs = self.obs_mp
        sigma = self.obs_sigma
        z = (model * obs)/sigma**2
        
        lp = sum(np.log(model) - 2*np.log(sigma) -  (model + obs)**2/(2*sigma**2) + np.log(i0(z)))
        return lp
