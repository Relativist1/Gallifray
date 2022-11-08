
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
import ehtim as eh
import sys
from gallifray.const import *


def ipole_to_obs(X, Y, Z, obs_ref, fov_m):
    ra = obs_ref.ra
    dec = obs_ref.dec
    rf = obs_ref.rf
    source = obs_ref.source 
    DX = fov_m
    DY = DX
    pix = len(X)
    NX= pix
    NY = NX
    Msun = 1.989e33
    MBH = M_sgr*Msun
    GNEWT = 6.6742e-8
    CL = 2.99792458e10
    L_unit = GNEWT * MBH  / (CL * CL)
    PC = 3.085678e18
    Dsource = D_sgr*PC
    psize = fov_m*L_unit/Dsource * 206264806247.1*1e-6*eh.RADPERAS/pix
    Zf = Z
    im1 = eh.image.Image((Zf),psize=psize,ra=ra, dec=dec, rf= rf, source=source)
    im1 = im1.rotate(np.pi/2)
    obs = im1.observe_same(obs_ref, add_th_noise=True,ttype='direct',verbose=False)
    blockp()

    obs.add_scans()
    obs = obs.avg_coherent(0.,scan_avg=True)
    vis_amp1 = obs.unpack(['amp'],conj=True)['amp']
    vis_err1 = obs.unpack(['sigma'],conj=True)['sigma'] 
    
    return vis_amp1, vis_err1

def grrt_to_obs(X, Y, Z, obs_ref, fov_m):
    ra = obs_ref.ra
    dec = obs_ref.dec
    rf = obs_ref.rf
    source = obs_ref.source 
    DX = fov_m
    DY = DX
    pix = len(X)
    NX= pix
    NY = NX
    Msun = 1.989e33
    MBH = M_sgr*Msun
    GNEWT = 6.6742e-8
    CL = 2.99792458e10
    L_unit = GNEWT * MBH  / (CL * CL)
    PC = 3.085678e18
    Dsource = D_sgr*PC
    psize = fov_m*L_unit/Dsource * 206264806247.1*1e-6*eh.RADPERAS/pix

    im1 = eh.image.Image(np.flipud(Z),psize=psize,ra=ra, dec=dec, rf= rf, source=source)

    obs = im1.observe_same(obs_ref, add_th_noise=True,ttype='direct',verbose=False)
    blockp()

    obs.add_scans()
    obs = obs.avg_coherent(0.,scan_avg=True)
    vis_amp1 = obs.unpack(['amp'],conj=True)['amp']
    vis_err1 = obs.unpack(['sigma'],conj=True)['sigma'] 
    
    return vis_amp1, vis_err1