#    Gaussian Model Class
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
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
"""
For Scattering Kernel
beam_params = [1.309*1e-6, 0.64*1e-6, 78*np.pi/180]
Reference :- Bower G.C. et. al., 2006, ApJ, 648, L127
"""


class sym_gauss(object):
    """Class for Symmetric Gaussian Model.
    
    Attributes:
        dim (int) : Dimenion along one axis (Square image)
        fov (int) : Field of view (in uas)
        pixel (int) : Pixel size
        x_off (float) : x offset from the center
        y_off (float) : y offset from the center
    """

    def __init__(self,I0,sigma, fov,dim=64,x_off=0, y_off=0):
        """Creates a Symmetric Gaussian Model.
        
        Args:
            I0 (float) : Total flux/intensity (Jy/pixel)
            sigma (float) : Width
            phi (float) : Orientation
            fov (int) : Field of view (in uas)
            dim (int) : Dimensions of the image along an axis (square image)
            x_off (float) : x offset from the center
            y_off (float) : y offset from the center
        Return:
            Disk model with parameters
        """
        
        self.I0 = I0
        self.sigma = sigma
        self.fov = fov
        self.dim = dim
        self.x_off = x_off
        self.y_off = y_off
        
        self.X = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.Y = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.psize = self.X[1] - self.Y[0]

    def sky_map(self):
        """Generates the intensity map of the model
            
        Returns:
            Intensity map of the model
        """
        psize = self.psize
        gauss_arr = np.zeros((self.dim, self.dim))
        for j, nj in enumerate(self.Y):
            for i, ni in enumerate(self.Y):
                gauss_arr[i-self.x_off][j-self.y_off] = np.exp(-(ni**2 + nj**2)/(2.*self.sigma**2))
        Gauss = gauss_arr*self.I0 #Normalise
        return Gauss
    
    def vis_data(self, fov, uv='default',A=-0.5,interp=None,points=512):
        """"Generate complex visibilites
            
        Return:
            vis_data: a dictionary object containing the information of complex visibilites and the baselines.
            
        """
        if uv == 'default':
            u1 = np.linspace(.1, self.fov, self.dim)
            v1 = np.linspace(.1, self.fov, self.dim)
            u = u1
            v = v1
            
        elif len(uv)==2:
            u = np.asarray(uv[0])
            v = np.asarray(uv[1])
            
        visibility = self.I0*np.exp(-(u**2- v**2)/(2*self.sigma))
        uv = np.sqrt(u**2 + v**2)
        
        bl_new = uv
        vis_n = visibility
    
        if interp:
            if not points:
                points = len(uv)
            bl_new = np.linspace(min(uv), max(uv), points)
            interp_vis = interp1d(np.asarray(uv), np.abs(visibility), kind=interp)
            vis3 = interp_vis(bl_new)
            vis_n = vis3


        vis_data = {'info': 'Complex Visibilites',
                  'vis' : vis_n,
                  'real': np.real(vis_n),
                  'imaginary': np.imag(vis_n),
                  'amp': np.abs(vis_n),
                  'u': u,
                  'v': v,
                  'bl': bl_new
                  }
        return vis_data
    
    def sky_blur(self, beam_params=[1.309, 0.64, 78*np.pi/180], beam_size=10):
        """Convolves the image with a gaussian kernel
           Default sets to Saggitarius A* Kernel
           Reference :- Bower G.C. et. al., 2006, ApJ, 648, L127
        Args:
            beam_params (list) : [fwhm_maj, fwhm_min, phi] where fwhm_maj (fwhm_min) is the full width half maximum of the major axis (minor axis) and phi is the position angle.
            beam_size (float) : Beam size
            
        Return:
            Blurred image convolved with the gaussian scattering kernel
        """
        if len(beam_params) != 3:
            raise Exception("beam_params must contain 3 values")
        
        image = self.sky_map()
        x,y = np.meshgrid(self.X,self.Y)
        
        fwhm_mj = beam_params[0]
        fwhm_min = beam_params[1]
        theta = beam_params[2]
        s_maj = fwhm_mj / (2. * np.sqrt(2. * np.log(2.)))
        s_min = beam_params[1] / (2. * np.sqrt(2. * np.log(2.)))
        x0 = np.cos(theta)
        y0 = np.sin(theta)
        
        Gauss = self.I0*np.exp(-(y*x0 + x*y0)**2/(2*(beam_size * s_maj)**2)- \
                          (x*x0 - y * y0)**2/(2.*(beam_size * s_min)**2))

        imarr_blur = fftconvolve(Gauss, image, mode='same')
        return imarr_blur

class asym_gauss(object):
    """Class for Asymmetric Gaussian Model.
    
    Attributes:
        dim (int) : Dimenion along one axis (Square image)
        fov (int) : Field of view (in uas)
        pixel (int) : Pixel size
        x_off (float) : x offset from the center
        y_off (float) : y offset from the center
    """

    def __init__(self,I0,sigma, A, phi, fov,dim=64,x_off=0, y_off=0):
        """Creates a Symmetric Gaussian Model.
        
        Args:
            I0 (float) : Total flux/intensity (Jy/pixel)
            sigma : Width
            A : Anisotropy parameter (0,1)
            phi (float) : Orientation
            fov (int) : Field of view (in uas)
            dim (int) : Dimensions of the image along an axis (square image)
            x_off (float) : x offset from the center
            y_off (float) : y offset from the center
        Return:
            Disk model with parameters
        """
        
        self.I0 = I0
        self.sigma = sigma
        self.s_min = self.sigma/np.sqrt(1+A)
        self.s_maj = self.sigma/np.sqrt(1-A)
        self.phi = phi
        self.fov = fov
        self.dim = dim
        self.x_off = x_off
        self.y_off = y_off
        
        self.X = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.Y = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.psize = self.X[1] - self.X[0]

    def sky_map(self):
        """Generates the intensity map of the model
                
            Returns:
                Intensity map of the model
        """

        x0 = np.cos(self.phi)
        y0 = np.sin(self.phi)
        gauss_arr = np.zeros((self.dim, self.dim))
        for j, nj in enumerate(self.Y):
            for i, ni in enumerate(self.Y):
                gauss_arr[i-self.x_off][j-self.y_off] = np.exp(-(nj*x0 + ni*y0)**2/(2*self.s_maj**2)-(ni*x0 - nj*y0)**2/(2.*self.s_min**2))
        Gauss = gauss_arr*self.I0 #Normalise
        return Gauss
    
    def vis_data(self, fov, uv='default', A=-0.5,interp=None,points=512):
        """"Generate complex visibilites
            
        Return:
            vis_data: a dictionary object containing the information of complex visibilites and the baselines.
            
        """
        if uv == 'default':
            u1 = np.linspace(.1, self.fov, self.dim)
            v1 = np.linspace(.1, self.fov, self.dim)
            u = u1*np.cos(self.phi) + v1*np.sin(self.phi)
            v = -u1*np.sin(self.phi) + v1*np.cos(self.phi)
            
        elif len(uv)==2:
            u = np.asarray(uv[0])
            v = np.asarray(uv[1])
            
        x0 = np.cos(self.phi)
        y0 = np.sin(self.phi)
        visibility = self.I0*np.exp(-(u*x0 + v*y0)**2/(2*self.s_maj**2)-(v*x0 - u*y0)**2/(2.*self.s_min**2))
        uv = np.sqrt(u**2 + v**2)
        
        bl_new = uv
        vis_n = visibility
        
        if interp:
            if not points:
                points = len(uv)
            bl_new = np.linspace(min(uv), max(uv), points)
            interp_vis = interp1d(np.asarray(uv), np.abs(visibility), kind=interp)
            vis3 = interp_vis(bl_new)
            vis_n = vis3

        vis_data = {'info': 'Complex Visibilites',
                  'vis' : vis_n,
                  'real': np.real(vis_n),
                  'imaginary': np.imag(vis_n),
                  'amp': np.abs(vis_n),
                  'u': u,
                  'v': v,
                  'bl': bl_new
                  }
        return vis_data
    
    def sky_blur(self, beam_params=[1.309, 0.64, 78*np.pi/180], beam_size=10):
        """Convolves the image with a gaussian kernel
           Default sets to Saggitarius A* Kernel
           Reference :- Bower G.C. et. al., 2006, ApJ, 648, L127
            Args:
                beam_params (list) : [fwhm_maj, fwhm_min, phi] where fwhm_maj (fwhm_min) is the full width half maximum of the major axis (minor axis) and phi is the position angle.
                beam_size (float) : Beam size
                
            Return:
                Blurred image convolved with the gaussian scattering kernel
        """
        
        if len(beam_params) != 3:
            raise Exception("beam_params must contain 3 values")
        
        image = self.sky_map()
        x,y = np.meshgrid(self.X,self.Y)
        
        fwhm_mj = beam_params[0]
        fwhm_min = beam_params[1]
        theta = beam_params[2]
        s_maj = fwhm_mj / (2. * np.sqrt(2. * np.log(2.)))
        s_min = beam_params[1] / (2. * np.sqrt(2. * np.log(2.)))
        x0 = np.cos(theta)
        y0 = np.sin(theta)
        
        Gauss = self.I0*np.exp(-(y*x0 + x*y0)**2/(2*(beam_size * s_maj)**2)- \
                          (x*x0 - y*y0)**2/(2.*(beam_size * s_min)**2))

        imarr_blur = fftconvolve(Gauss, image, mode='same')
        return imarr_blur
