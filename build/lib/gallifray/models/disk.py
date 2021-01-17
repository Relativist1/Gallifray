from __future__ import division
from __future__ import print_function

from builtins import object
import numpy as np
from scipy.signal import fftconvolve
from scipy.interpolate import interp1d
from scipy.special import jv
"""
For Scattering Kernel
beam_params = [1.309*1e-6, 0.64*1e-6, 78*np.pi/180]
Reference :- Bower G.C. et. al., 2006, ApJ, 648, L127
"""


class disk(object):
    """Class for Disk Model.
    
    Attributes:
        dim (int) : Dimenion along one axis (Square image)
        fov (float) : Field of view (in uas)
        pixel (int) : Pixel size
        x_off (float) : offset in x
        y_off (float) : offset in y
    """

    def __init__(self,I0, R, fov,dim=512,x_off=0, y_off=0):
        """Creates a Disk Model.
        
        Args:
            I0 (float) : Total flux/intensity (Jy/pixel)
            R (float) : Outer radius (in uas)
            fov (int) : Field of view (in uas)
            dim (int) : Dimensions of the image along an axis (square image)
            x_off (float) : x offset from the center
            y_off (float) : y offset from the center
        Return:
            Disk model with parameters
        """
        
        self.I0 = I0
        self.R = R
        self.fov = fov
        self.dim = dim
        self.x_off = x_off
        self.y_off = y_off
        
        self.X = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.Y = np.linspace(-self.fov/2, self.fov/2,self.dim)
        self.pixel = self.X[1] - self.X[0]

    def sky_map(self):
        """Generates the intensity map of the model

        Returns:
            Intensity map of the model
        """
        disk = np.zeros((self.dim, self.dim))
        for j, nj in enumerate(self.Y):
            for i, ni in enumerate(self.Y):
                R1 = nj**2 + ni**2
                if R1 < self.R**2:
                    disk[i-self.x_off][j-self.y_off] = self.I0
        return disk
    
    def vis_data(self, fov, uv='default',A=-0.5,interp=None,points=512):
        """"Generate complex visibilites
            
        Return:
            vis_data: a dictionary object containing the information of complex visibilites and the baselines.
            
        """
        if uv == 'default':
            u = np.linspace(.1, self.fov, self.dim)
            v = np.linspace(.1, self.fov, self.dim)
        elif len(uv)==2:
            u = np.asarray(uv[0])
            v = np.asarray(uv[1])
            
        rho = np.sqrt(u**2 + v**2)
        visibility = self.R*jv(1,(2*np.pi*rho)*self.R)/rho
        uv = np.sqrt(u**2 + v**2)
        
        if interp!=None and interp=='spline':
            def cubic_spline_interp(x,y,new_x,a=-0.5) :
                delta = x[1]-x[0]
                F = np.zeros(len(new_x))
                for j in range(len(new_x)) :
                    dx = (x-new_x[j])/delta
                    weight = (-a)*(dx**3+5*dx**2+8*dx+4)*(dx>=-2)*(dx<-1) + \
                        (-(a+2)*dx**3-(a+3)*dx**2+1)*(dx>=-1)*(dx<=0) + \
                        ((a+2)*dx**3-(a+3)*dx**2+1)*(dx>0)*(dx<=1) + \
                        (-a)*(-dx**3+5*dx**2-8*dx+4)*(dx>1)*(dx<=2)
                    F[j] = np.sum(weight*y)
                return F
            
            bl_new = np.linspace(min(uv), max(uv), points)
            vis3 = cubic_spline_interp(uv,np.abs(visibility),bl_new,a=A)
            vis_n = vis3/max(vis3)

        if interp!=None and interp!='spline':
            if not points:
                points = len(uv)
            bl_new = np.linspace(min(uv), max(uv), points)
            interp_vis = interp1d(np.asarray(uv), np.abs(visibility), kind=interp)
            vis3 = interp_vis(bl_new)
            vis_n = vis3/max(vis3)

        if interp==None:
            bl_new = uv
            vis_n = visibility/max(visibility)


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

