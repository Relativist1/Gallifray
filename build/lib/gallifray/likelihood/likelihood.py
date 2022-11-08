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
from gallifray.likelihood.image_obs import *
from gallifray.const import *
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
        obs_data.add_scans()
        obs_data = obs_data.avg_coherent(0.0,scan_avg=True)
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
        if self.model[1]=='sym_gauss':
            I0, sigma = self.param
            dim = len(self.obs_amp)
            imarr = sym_gauss(I0, sigma, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='asym_gauss':
            I0, A, sigma, phi = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='disk':
            I0, R = self.param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='crescent':
            I0, R, psi, tau, phi = self.param
            dim = len(self.obs_amp)
            imarr = crescent(I0, R, psi, tau, phi,fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='xsring':
            I0, R_p, R_n, ecn, f, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=64)
            IM = imarr.sky_map()

        ra = self.obs_data.ra
        dec = self.obs_data.dec
        rf = self.obs_data.rf
        source = self.obs_data.source
        # fov_m = fov
        # DX = fov
        # DY = DX
        # pix = len(self.obs_amp)
        # NX= pix
        # NY = NX
        # Msun = 1.989e33
        # MBH = M_sgr*Msun
        # GNEWT = 6.6742e-8
        # CL = 2.99792458e10
        # L_unit = GNEWT * MBH  / (CL * CL)
        # PC = 3.085678e18
        # Dsource = D_sgr*PC
        # psize = fov_m*L_unit/Dsource * 206264806247.1*1e-6*eh.RADPERAS/pix
        im1 = eh.image.Image(IM/np.sum(IM)*I0,psize=imarr.psize*eh.RADPERUAS,ra=ra, dec=dec, rf= rf, source=source)
        obs1 = im1.observe_same(self.obs_data, add_th_noise=True,ttype='direct',verbose=False)
        blockp()

        obs1.add_scans()
        obs1 = obs1.avg_coherent(0.,scan_avg=True)
        model_amp = obs1.unpack(['amp'],conj=True)['amp']

        model = model_amp/np.max(model_amp)
        obs = self.obs_amp/np.max(self.obs_amp)
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
        if self.model[1]=='sym_gauss':
            I0, sigma = self.param
            dim = len(self.obs_amp)
            imarr = sym_gauss(I0, sigma, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='asym_gauss':
            I0, A, sigma, phi = self.param
            dim = len(self.obs_amp)
            imarr = asym_gauss(I0, A, sigma, phi, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='disk':
            I0, R = self.param
            dim = len(self.obs_amp)
            imarr = disk(I0, R, fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='crescent':
            I0, R, psi, tau, phi = self.param
            dim = len(self.obs_amp)
            imarr = crescent(I0, R, psi, tau, phi,fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='xsring':
            I0, R_p, R_n, ecn, f, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsring(I0, R_p, R_n, ecn, f, phi,fov, dim=64)
            IM = imarr.sky_map()
            
        if self.model[1]=='xsringauss':
            I0, R_p, R_n, ecn, f, gax, aq, gq, phi = self.param
            dim = len(self.obs_amp)
            imarr = xsringauss(I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov, dim=64)
            IM = imarr.sky_map()


        ra = self.obs_data.ra
        dec = self.obs_data.dec
        rf = self.obs_data.rf
        source = self.obs_data.source
        # fov_m = fov
        # DX = fov
        # DY = DX
        # pix = len(self.obs_amp)
        # NX= pix
        # NY = NX
        # Msun = 1.989e33
        # MBH = M_sgr*Msun
        # GNEWT = 6.6742e-8
        # CL = 2.99792458e10
        # L_unit = GNEWT * MBH  / (CL * CL)
        # PC = 3.085678e18
        # Dsource = D_sgr*PC
        # psize = fov_m*L_unit/Dsource * 206264806247.1*1e-6*eh.RADPERAS/pix
        im1 = eh.image.Image(IM/np.sum(IM)*I0,psize=imarr.psize*eh.RADPERUAS,ra=ra, dec=dec, rf= rf, source=source)
        obs1 = im1.observe_same(self.obs_data, add_th_noise=True,ttype='direct',verbose=False)
        blockp()

        obs1.add_scans()
        obs1 = obs1.avg_coherent(0.,scan_avg=True)
        model_amp = obs1.unpack(['amp'],conj=True)['amp']

        model = model_amp/np.max(model_amp)
        obs = self.obs_amp/np.max(self.obs_amp)
        sigma = self.obs_sigma
        z = (model * obs)/sigma**2
        lp = sum(np.log(model) - 2*np.log(sigma) -  (model + obs)**2/(2*sigma**2) + np.log(i0(z)))/(dim + len(self.param))
        return lp

    def ln_gaussian_physical(self, pr_param, exec_c, fov_m, code_type=None):
        """GALLIFRAY :: Gaussian distribution likelihood for general purposes
        
        Args:
            exec_c (list): executable along with extra arguments
        Return:
            Calculated Likelihood Distribution
        """

        X, Y, Z = prep_grrt(self.param, exec_c, fov_m) 
        if code_type=='ipole':
            model_amp, err = ipole_to_obs(X, Y, Z, self.obs_data, fov_m)
        else:
            model_amp, err = grrt_to_obs(X, Y, Z, self.obs_data, fov_m)
            
        model = model_amp/np.max(model_amp)
        obs = self.obs_amp/np.max(self.obs_amp)
        sigma = self.obs_sigma
        
        lp = -0.5*sum((model - obs)**2/sigma + np.log(sigma) + np.log(2*np.pi))


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