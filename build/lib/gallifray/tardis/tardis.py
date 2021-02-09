#    Tardis Class - Triangular Distribution Plotting for MCMC sampling analysis
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
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 18


def Tardis(samples,
                  labels=None,
                  truths=None,
                  savefig=None,
                  contour_levels=5,
                  shade=True,
                  diag_shade=True,
                  shade_color='Darkblue',
                  diag_shade_color = 'Darkblue',
                  truth1d=True,
                  truth2d=False,
                  truth_titles=False,
                  color_truth='k',
                  lw_truth=1.5,
                  lw_1d=2,
                  fontsize=20,
                  pad_inches=0.1,
                  dpi=500,
                  **kwargs):
    """Tardis: Triangle-corner distribution plotting for MCMC sampling analysis.C
    """
    
    dim = len(samples.T)
    whspace = 0.15
    plotdim = 3*(dim + (dim - 1.0)*whspace)
    size = 1.5 + plotdim + 0.6
    
    fig, axes = plt.subplots(dim,dim,figsize=(size,size))
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xticks(size=14)
    plt.yticks(size=14)
    
    for i in range(dim):
        for j in range(i):
            ax = axes[i, j]
            ax.xaxis.set_major_locator(plt.MaxNLocator(4))
            ax.yaxis.set_major_locator(plt.MaxNLocator(4))
            fig.delaxes(axes[j][i])
            ax.tick_params(labelbottom=False, labelleft=False, labelright=False,
                           bottom=True, top=False, left=True, right=False,
                           direction='in', labelsize='large')
            if (j==0) :
                ax.tick_params(labelleft=True)
                if labels is not None:
                    ax.set_ylabel(labels[i],fontsize=fontsize)
            if (i==dim-1) :
                ax.tick_params(labelbottom=True)
                if labels is not None:
                    ax.set_xlabel(labels[j],fontsize=fontsize)
            if (j==dim-1) :
                ax.tick_params(labelright=True)
            if truth2d is not False:
                ax.axvline(truths[i],color=color_truth,lw=lw_truth)
                ax.axhline(truths[i],color=color_truth,lw=lw_truth)
                
            M = sns.kdeplot(x=samples[:,j], y=samples[:,i],ax=ax, color=shade_color, shade=shade, levels=contour_levels,
                            **kwargs)

    for i in range(dim):
        ax = axes[i, i]
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax.tick_params(labelleft=False, labelright=False, labelbottom=False,
                       labeltop=False, bottom=True, top=True, left=True,
                       right=True, direction='in', length=4, labelsize='large')
        
        if i==0:
            if labels is not None:
                ax.set_xlabel(labels[-1],fontsize=fontsize)
                ax.set_ylabel(labels[i],fontsize=fontsize)
        if i==dim-1:
            ax.tick_params(labelbottom=True)
            if labels is not None:
                ax.set_xlabel(labels[-1],fontsize=fontsize)
        
        if truth1d is not False:
            ax.axvline(truths[i],color=color_truth,lw=lw_truth)
        
        N = sns.kdeplot(x=samples[:,i],ax=ax,shade=diag_shade, color=diag_shade_color, lw=lw_1d, **kwargs)
            
    if savefig is not None:
        fname = savefig
        plt.savefig(fname, bbox_inches='tight', pad_inches = pad_inches,dpi=dpi)
        
    return M, N
    
