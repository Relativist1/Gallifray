#    Priors Class
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

import numpy as np

models_list = ["sym_gaussian", "asym_gaussian", "disk", "crescent", "xsring", "xsringauss"]

pi = 3.141592653589793


class priors():
    """Prior Class
        
    """
    def __init__(self, param):
        """
        
        Args:
            param (list): list of parameters
            
        Return:
            
        """
        self.param = param
    
    def lnprior(self, model_type, type='uniform'):
        """
        Args:
            model_type (str): model
            type: type of prior
                
                'flat' : Set uniform priors
                'gaussian' : Set gaussian priors
            
        Return:
            
        """
        if model_type not in models_list:
            raise("model not yet implemented!")
    
        if model_type=='sym_gaussian':
            pr = prior_sym_gauss(self.param)
        if model_type=='asym_gaussian':
            pr = prior_asym_gauss(self.param)
        if model_type=='disk':
            pr = prior_disk(self.param)
        if model_type=='crescent':
            pr = prior_crescent(self.param)
        if model_type=='xsring':
            pr =  prior_xsring(self.param)
        if model_type=='xsringauss':
            pr = prior_xsringauss(self.param)
            
        return pr
    
    
def prior_sym_gauss(param):
       
    I0, S, phi = param

    if 0<I0<10 and 0<S<100 and 0<phi<np.pi:
        return 0
    else:
        return -np.inf

def prior_asym_gauss(param):
    
    I0, S, A, phi = param
    
    if 0<I0<5 and 0<S<100 and 0<A<0.999 and 0<phi<np.pi:
        return 0
    else:
        return -np.inf

def prior_disk(param):
    
    I0, R_p, = param
    if 0<I0<100 and 0<R_p<100:
        return 0
    else:
        return -np.inf

def prior_crescent(param):
    
    I0, R_p, R_n, ecn, phi = param
    if 0<I0<100 and 0<R_p<100 and 0<R_n<100 and -1<ecn<1 and -np.pi<phi<np.pi:
        return 0
    else:
        return -np.inf


def prior_xsring(param):
    
    I0, R_p, R_n, ecn, f, phi = param
    if 0<I0<100 and 0<R_p<100 and 0<R_n<100 and -1<ecn<1 and 0<f<1 and -np.pi<phi<np.pi:
        return 0
    else:
        return -np.inf
        
def prior_xsringauss(param):

    I0, R_p, R_n, ecn, f, gax, aq, gq, phi = param
    if 0<I0<5 and 0<R_p<100 and 0<R_n<100 and -1<ecn<1 and 0<f<1 and 0<gax<1 and 0<aq<1 and 0<gq<1 and -np.pi<phi<np.pi:
        return 0
    else:
        return -np.inf


#def lnprior(theta):
#    a, b, c = theta
#    #flat priors on b, c
#    if not 1.0 < b < 2.0 and 1.0 < c < 2.0:
#        return -np.inf
#    #gaussian prior on a
#    mu = 10
#    sigma = 1
#    return np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(a-mu)**2/sigma**2
