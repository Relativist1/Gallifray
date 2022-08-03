"""
.. module:: gallifray
    :synopsis: Geometric Modelling and Parameter Estimation framework

.. moduleauthor:: Saurabh (sbhkmr1999@gmail.com)

"""
from __future__ import division
from __future__ import print_function



models_list = ['sym_gauss', 'asym_gauss', 'disk', 'crescent', 'xsring', 'xsringauss']


label_sym_gauss =  [r"$I_0$", r"$\sigma(\mu as)$", r"${\phi}$"]
label_asym_gauss =  [r"$I_0$", r"$\sigma(\mu as)$", "$A$", r"${\phi}$"]
label_disk =        [r"$I_0$", r"$R(\mu as)$", r"${\phi}$"]
label_crescent =    [r"$I_0$", r"$R_p (\mu as)$", "$\psi$", r"${\tau}$",  r"${\phi}$"]
label_xsring =      [r"$I_0$", r"$R_p (\mu as)$", "$R_n (\mu as)$", r"${\epsilon}$", r"f", r"${\phi}$"]
label_xsringauss =  [r"$I_0$", r"$R_p (\mu as)$", "$R_p (\mu as)$", r"${\epsilon}$", r"f",r"${g_{ax}}$",r"${a_{q}}$",r"${g_{q}}$", r"${\phi}$"]

def Labels(model):
    if model not in models_list:
        raise Exception("Supplied model error, check again!")
    if model==models_list[0]:
        return label_sym_gauss
    if model==models_list[1]:
        return label_asym_gauss
    if model==models_list[2]:
        return label_disk
    if model==models_list[3]:
        return label_crescent
    if model==models_list[4]:
        return label_xsring
    if model==models_list[5]:
       return label_xsringauss
