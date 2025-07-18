import batman, catwoman
import os, sys
import pickle
from datetime import datetime
import numpy as np
import matplotlib
from matplotlib import cm
from matplotlib.patches import Arrow
from matplotlib.gridspec import GridSpec
from astropy.table import Table, Column, vstack
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns 
sns.set_style('ticks')
import juliet

def load_plt_params():
    """ Load in plt.rcParams and set (based on paper defaults).
    """
    params = Table.read('rcParams.txt', format='csv')
    for i, name in enumerate(params['name']):
        try:
            plt.rcParams[name] = float(params['value'][i])
        except:
            plt.rcParams[name] = params['value'][i]
    return params

pltparams = load_plt_params()

fig = plt.figure(figsize=(25,24))
fig.set_facecolor('w')

# First, cut the parts where the spectra will go:
spacing = 15
gs = GridSpec(400+spacing*3, 1, figure=fig, wspace = 0.1)
ax = {}

for name, index in zip(['003', '001', '104', '002'], [0, 1, 2, 3]):

    if index == 0:

        ax[name] = fig.add_subplot(gs[:100, 0])

    else:

        ax[name] = fig.add_subplot(gs[index*(100+spacing):index*(100+spacing) + 100, 0])

    # Set y-labels:
    ax[name].set_ylabel('Relative Depth (ppm)')

# Set x-labels:
for name in ['003', '001', '104', '002']:

    # Take the chance to remove ticks from top plots in the bottom:
    if name != '002':

        ax[name].tick_params(labelbottom=False)

    else:

        ax[name].set_xlabel(r'Wavelength ($\mu$m)')

# Set limits of plots:
for name in ['003', '001', '104', '002']:

    ax[name].set_xlim(0.6, 5.0)
    ax[name].set_ylim(-600.,600.)

# All right --- now we are ready to plot things. Let's do a for loop to generate the plots:
espinoza_color = 'black'#'cornflowerblue'
gressier_color = '#3d3b8e'#'darkblue'
grant_color = '#4e937a'#'forestgreen'
stevenson_color = '#f0544f'#'firebrick'
canas_color = '#d81e5b'#'orangered'
ms = 13
ew = 2
alpha = 0.8
for name, date in zip(['003', '001', '104', '002'], ['June 22', 'June 28', 'July 23', 'October 28']):

    ax[name].text(4.2, 400,date+', 2023', fontsize = 30)

    w_ne, d_ne, derr_ne = np.loadtxt('spectra/JWST_PRISM_Espinoza'+name+'.txt', unpack = True, usecols = (0, 3, 4))
    w_ag, d_ag, derr_ag = np.loadtxt('spectra/JWST_PRISM_Gressier'+name+'.txt', unpack = True, usecols = (0, 3, 4))
    w_dg, d_dg, derr_dg = np.loadtxt('spectra/JWST_PRISM_Grant'+name+'.txt', unpack = True, usecols = (0, 3, 4))
    w_scks, d_scks, derr_scks = np.loadtxt('spectra/JWST_PRISM_Stevenson'+name+'.txt', unpack = True, usecols = (0, 3, 4))
    w_cc, d_cc, derr_cc = np.loadtxt('spectra/JWST_PRISM_Canas'+name+'.txt', unpack = True, usecols = (0, 3, 4))

    ax[name].errorbar(w_ne, d_ne - np.nanmedian(d_ne), derr_ne, fmt = 'o', \
                      ecolor = espinoza_color, \
                      mfc = espinoza_color, \
                      mec = espinoza_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'NE', alpha = alpha) 

    ax[name].errorbar(w_ag, d_ag - np.nanmedian(d_ag), derr_ag, fmt = 'o', \
                      ecolor = gressier_color, \
                      mfc = gressier_color, \
                      mec = gressier_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'AG', alpha = alpha)

    ax[name].errorbar(w_dg, d_dg - np.nanmedian(d_dg), derr_dg, fmt = 'o', \
                      ecolor = grant_color, \
                      mfc = grant_color, \
                      mec = grant_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'DG', alpha = alpha)

    ax[name].errorbar(w_scks, d_scks - np.nanmedian(d_scks), derr_scks, fmt = 'o', \
                      ecolor = stevenson_color, \
                      mfc = stevenson_color, \
                      mec = stevenson_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'SC/KS', alpha = alpha)

    ax[name].errorbar(w_cc, d_cc - np.nanmedian(d_cc), derr_cc, fmt = 'o', \
                      ecolor = canas_color, \
                      mfc = canas_color, \
                      mec = canas_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'CC', alpha = alpha)

    #ax[name].set_xscale('log')

ax['003'].legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc='lower left', \
           mode='expand', borderaxespad=0, ncol=5, fontsize = 30, frameon=False)

plt.tight_layout()
plt.savefig('tspectra_all.pdf',
            dpi=250,#, rasterize=True,
            #transparent=True,
            bbox_inches='tight')
