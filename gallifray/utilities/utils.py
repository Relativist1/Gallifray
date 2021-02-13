#    Utitlity functions
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
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import griddata
from scipy.fftpack import fft2, fftshift

import sys, os
from contextlib import contextmanager

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout

try:
    with suppress_stdout():
        import ehtim as eh
except ImportError:
    print("ehtim not installed!")


def save_fits(im, fname, psize):
    """Save image data to a fits file.

       Args:
            im (array): 2-d image data 
            fname (str): Image file name for saving

       Returns:
            saved fits file
    """
    mjd = 51544
    #Set for M87
    ra = 12+ 30/60+ 49/3600
    dec = 12 + 23/60 + 28/3600
    
    dim = im.shape[0]
    deg = np.pi/180
    hr = 15.0*deg
    radp = deg/3600.0
    #psize = 0.000325*radp/1024
    rf = 227000000000 #freq
    
    header = fits.Header()
    header['OBJECT'] = 'M87'
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CDELT1'] = -psize/deg
    header['CDELT2'] = psize/deg
    header['OBSRA'] = ra * 180/12.
    header['OBSDEC'] = dec
    header['FREQ'] = rf 
    header['CRPIX1'] = dim/2 + 0.5
    header['CRPIX2'] = dim/2 + 0.5

    image = np.reshape(im, (dim, dim))[::-1, :]
    hdu = fits.PrimaryHDU(image, header=header)
    hdulist = fits.HDUList([hdu])
    
    hdulist.writeto(fname, overwrite=True)

    return
    
def save_vamp_fits(im, fname):
    """Computes visibility amplitude of the imagea and saves as fits files.

       Args:
            im (array): 2-d image data
            fname (str): Image file name for saving

       Returns:
            saved fits file
    """
    img = fits.open(im)
    img_data = img[0].data
    
    mjd = 51544
    #Set for M87
    ra = 12+ 30/60+ 49/3600
    dec = 12 + 23/60 + 28/3600

    dim = im.shape[0]
    deg = np.pi/180
    hr = 15.0*deg
    radp = deg/3600.0
    psize = 0.000325*radp/1024
    rf = 227000000000 #freq

    header = fits.Header()
    header['OBJECT'] = 'M87'
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CDELT1'] = -psize/DEGREE
    header['CDELT2'] = psize/DEGREE
    header['OBSRA'] = ra * 180/12.
    header['OBSDEC'] = dec
    header['FREQ'] = rf
    header['CRPIX1'] = dim/2 + 0.5
    header['CRPIX2'] = dim/2 + 0.5
    
    vis_amp = np.abs(fftshift(fft2(img_data))).flatten()
    hdu = fits.PrimaryHDU(vis_amp, header=header)
    hdulist = fits.HDUList([hdu])

    hdulist.writeto(fname, overwrite=True)

    return

def raytraced_image(fname,show=False):
    """Return the raytraced image.

       Args:
            fname (str): Image file file
            show (bool): (optional) If true, shows the plot

       Returns:
            x_coordinate: 2-d array
            y_coordinate: 2-d array
            image: 2-d array
    """
    data = np.loadtxt(fname,unpack=True)
    xcoord = data[4]
    ycoord = data[5]
    inten = data[6]

    #Interpolation
    res  = 1024
    new_x = np.linspace(-30,30,res)
    new_y = np.linspace(-30,30,res)
    mesh_x,mesh_y = np.meshgrid(new_x,new_y)
    image = griddata((xcoord,ycoord),inten,(mesh_x,mesh_y),method='nearest')
    if show:
        plt.pcolormesh(mesh_x,mesh_y,image,cmap='afmhot',shading='gouraud')

    image_data = [mesh_x , mesh_y, image]
    return image

def obs_uv(im_fits, out_file,array='EHT2017', t_int=5, t_adv=600, t_start=0, t_stop=24, bw=4e9, thermal_noise=True, save=False):
    """Returns the observation file over the EHT baselines (8 working stations).
    
        Args:
             im_fits (str): Image file file
             out_file (str): Output file name
             t_int (int):  Integration time in seconds (seconds)
             t_adv (int): advance time between scans (seconds)
             t_start (float): GMST time of the start of the observation (hours)
             t_stop (float): GMST time of the end of the observation (hours)
             bw (float): bandwidth in Hz
             thermal_noise (bool): if True, add thermal noise to the observed data
        Returns:
             Observed uvfits file
    """
    fname_arr = '/Users/geodesix/eht-imaging/arrays/EHT2017.txt'
    
    im = eh.image.load_image(im_fits)
    arr = eh.array.load_txt(fname_arr)
    if array!='EHT2017':
        arr = eh.array.load_txt(array)
    
    if not thermal_noise:
        obs = im.observe(arr, t_int, t_adv, t_start, t_stop, bw,
                        sgrscat=False, ampcal=True, phasecal=False,
                        ttype='fast', add_th_noise=False)
    else:
        obs = im.observe(arr, t_int, t_adv, t_start, t_stop, bw,
                        sgrscat=False, ampcal=True, phasecal=False,
                        ttype='fast', add_th_noise=True)
    if not save:
        obs.save_uvfits(out_file)
    return obs
