#    Xsringauss Model Class
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
from scipy.special import jv
from scipy.interpolate import interp1d
"""
For Scattering Kernel
beam_params = [1.309*1e-6, 0.64*1e-6, 78*np.pi/180]
Reference :- Bower G.C. et. al., 2006, ApJ, 648, L127
"""

class xsringauss(object):
    """Class for xsringauss Model.
    Reference: Benkevitch, L., Akiyama, K., Lu, R., Doeleman, S., & Fish, V. 2016,
    arXiv:1609.00055
    
    Attributes:
        dim (int) : Dimenion along one axis (Square image)
        fov (float) : Field of view (in uas)
        pixel (int) : Pixel size
        x_off (float) : offset in x
        y_off (float) : offset in y
    """

    def __init__(self,I0, R_p, R_n, ecn, f, gax, aq, gq, phi,fov,dim=512,x_off=0, y_off=0):
        """Creates a xsringauss Model.
        
        Args:
            I0 (float) : Total flux/intensity (Jy/pixel)
            R_p (float) : Outer radius (in uas)
            R_n (float) : Inner radius (in uas)
            ecn (float) : Eccentricity (0,1)
            f (float) : Fading parameter (0,1)
            gax (float) : FWHM of main axis
            aq (float) : Axial ratio
            gq (float) : Fraction of gaussian flux in total flux (0,1)
            phi (float) : Orientation
            fov (int) : Field of view (in uas)
            dim (int) : Dimensions of the image along an axis (square image)
            x_off (float) : x offset from the center
            y_off (float) : y offset from the center

        Return:
            xsringauss model with parameters
        """
        
        self.I0 = I0
        self.R_p = R_p
        self.R_n = R_n
        self.ecn = ecn
        self.f = f
        self.gax = gax
        self.aq =  aq
        self.gq = gq
        self.phi = phi
        self.fov = fov
        self.dim = dim
        self.x_off = x_off
        self.y_off = y_off
        
        self.d = self.ecn * (self.R_p - self.R_n)
        self.a = self.gax * self.R_p
        self.b = self.aq*self.a
        
        self.X = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.Y = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.psize = (self.X[1]-self.X[0])
        
    def sky_map(self):
        """Generates the intensity map of the model
         
        Returns:
            Intensity map of the model
        """
    
        xsringauss_arr = np.zeros((self.dim, self.dim))
        for j, nj in enumerate(self.Y):
            for i, ni in enumerate(self.Y):
                R1 = (nj - self.d)**2 + (ni-self.d)**2
                R2 = nj**2 + ni**2
                x0 = ni*np.cos(self.phi) - nj*np.sin(self.phi)
                y0 = nj*np.cos(self.phi) + ni*np.sin(self.phi)
                if R1 > self.R_n**2 and R2 < self.R_p**2:
                    xs0 = (2*self.I0*(1-self.gq)/np.pi)/(self.R_p**2-self.R_n**2*(1+self.d/self.R_p) \
                                                         - (1-self.f)*(self.d*self.R_n**2/self.R_p))
                    xsring = xs0*((1/2)*self.f*(1.0-x0/self.R_p) + (1/2)*(1.0+x0/self.R_p))
                    xsgauss0 = self.gq*self.I0/(2.*np.pi*self.a*self.b)
                    Igauss = xsgauss0*np.exp(-(1/2)*((x0-(self.d-self.R_n))/self.a)**2 -0.5*((y0/self.b)**2))
                    xsringauss_arr[i-self.x_off][j-self.y_off] = xsring + Igauss

        return xsringauss_arr
    
    def vis_data(self, fov, uv='default',A=-0.5,interp=None,points=512):
        """"Generate complex visibilites
            
        Return:
            vis_data: a dictionary object containing the information of complex visibilites and the baselines.
            
        """
        if uv == 'default':
            u1 = np.linspace(.1, self.fov, self.dim)
            v1 = np.linspace(.1, self.fov, self.dim)
        elif len(uv)==2:
            u1 = np.asarray(uv[0])
            v1 = np.asarray(uv[1])
            
        u = u1*np.cos(self.phi) + v1*np.sin(self.phi)
        v = -u1*np.sin(self.phi) + v1*np.cos(self.phi)
        rho = np.sqrt(u**2 + v**2)
        d = self.ecn*(self.R_p - self.R_n)
        h = (2/np.pi)*(1/((self.R_p**2 - self.R_n**2*(1+self.d/self.R_p))*self.f + \
                          (self.R_p**2 - self.R_n**2*(1-self.d/self.R_p))))
        h0 = self.f*h
        xi =  (h + h0)/2
        eta = (h - h0)/2
        circ1 = self.R_p*jv(1,(2*np.pi*rho)*self.R_p)
        circ2 = self.R_n*jv(1,(2*np.pi*rho)*self.R_n)
        
        dcirc1 = u*self.R_p*(np.pi*self.R_p*(jv(0,(2*np.pi*rho)*self.R_p) - \
                                             jv(2,(2*np.pi*rho)*self.R_p))/rho**2 - \
                           jv(1,(2*np.pi*rho)*self.R_p)/rho**3)
        dcirc2 = u*self.R_n*(np.pi*self.R_n*(jv(0,(2*np.pi*rho)*self.R_n) - \
                                             jv(2,(2*np.pi*rho)*self.R_n))/rho**2 - \
                           jv(1,(2*np.pi*rho)*self.R_n)/rho**3)
        
        vis1 = xi*circ1 + eta*(1j/2*np.pi)*dcirc1*u
        vis2 = xi - eta*(self.d/self.R_p)*circ1 + eta*(self.R_n/self.R_p)*(1j/2*np.pi)*dcirc2*u
        rho2= (u*self.a)**2 + (v*self.b)**2
        k = 1/(2*np.sqrt(2*np.log(2)))
        v_r = vis1 - np.exp(1j*2*np.pi*d*u*vis2)
        v_g = np.exp(-2*(np.pi*k)**2*rho2 - 1j*2*np.pi*u*(self.R_n-self.d))
        visibility = self.I0*((1-self.gq)*v_r + self.gq*v_g)
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

        Gauss = self.I0*np.exp(-(y * x0 + x * y0)**2/(2*(beam_size * s_maj)**2)- \
                          (x * x0 - y * y0)**2/(2.*(beam_size * s_min)**2))

        imarr_blur = fftconvolve(Gauss, image, mode='same')
        return imarr_blur

