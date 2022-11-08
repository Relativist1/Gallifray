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

models_list = ["sym_gauss", "asym_gauss", "disk", "crescent", "xsring", "xsringauss", "physical"]

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
    
    def lnprior(self, model_type, p_type='uniform', use_default=True):
        """
        Args:
            model_type (str): model (eg. ['geom', 'sym_gauss'] or ['geom', 'xsring'] or ['physical',model])
            p_type: type of prior
                
                'uniform' : Set uniform priors
                'gaussian' : Set gaussian priors
            
        Note:
            For physical model extra args for parameters priors (type:dict) : {str : [mean, min, max, sigma]}

            eg: pr_params = {'a'   : [0, -1, 1, 0.1],
                            'inc'  : [0, -90, 90, 0.1],
                             'w'   : [0, -5, 5, 0.1]}
            
        """

        # if isinstance(model_type[1],str)==True:
        #     if model_type[1] not in models_list:
        #         raise("model not yet implemented!")
        

        if model_type[0]=='geom':
            
            if not use_default:
                pr = prior_gen(p0, param, p_type='uniform')
            else:
                if model_type[1]=='sym_gauss':
                    pr = prior_sym_gauss(self.param)
                if model_type[1]=='asym_gauss':
                    pr = prior_asym_gauss(self.param)
                if model_type[1]=='disk':
                    pr = prior_disk(self.param)
                if model_type[1]=='crescent':
                    pr = prior_crescent(self.param)
                if model_type[1]=='xsring':
                    pr =  prior_xsring(self.param)
                if model_type[1]=='xsringauss':
                    pr = prior_xsringauss(self.param)
            
        elif model_type[0]=='physical':
            pr = prior_gen(self.param, model_type[1],p_type)
        return pr 

def ln_g(x, mu, sigma):
    return np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(x-mu)**2/sigma**2

def prior_sym_gauss(param):
       
    I0, S = param

    if 0<I0<10 and 0<S<100:
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

def prior_gen(p0, param, p_type='uniform'):
    
    m1 = []
    m2 = []
    out = []
    
    mu = []
    sigma = []

    # for key, value in dict(pr_params).items():
    #     if value == 'freeze':
    #         pr_params.pop(key)
            
    #stores the min/max values

    for i in param.keys():
        mu.append(param[i][0])
        m1.append(param[i][1])
        m2.append(param[i][2])
        sigma.append(param[i][3])

    #checks the params if they are within limits
            
    for i in range(len(p0)):
        if m1[i]<p0[i]<m2[i]:
            out.append(0)

        else:
            out.append(-np.inf)
    #checks if any inf value
    final = 0
    if p_type=='gaussian':
        for i in range(len(p0)):
            final += ln_g(p0[i],mu[i],sigma[i])

    if min(out)==-np.inf:
        return -np.inf
    else:
        return final

