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

fig = plt.figure(figsize=(15,24))
fig.set_facecolor('w')

# First, cut the parts where the spectra will go:
spacing = 15
gs = GridSpec(400+spacing*3, 1, figure=fig, wspace = 0.1)
ax = {}

for name, index in zip(['visit1', 'visit2'], [0, 1]):

    if index == 0:

        ax[name] = fig.add_subplot(gs[:100, 0])

    else:

        ax[name] = fig.add_subplot(gs[index*(100+spacing):index*(100+spacing) + 100, 0])

    # Set y-labels:
    ax[name].set_ylabel('Transit depth (ppm)')

# Set x-labels:
for name in ['visit1', 'visit2']:

    # Take the chance to remove ticks from top plots in the bottom:
    if name != 'visit2':

        ax[name].tick_params(labelbottom=False)

    else:

        ax[name].set_xlabel(r'Wavelength ($\mu$m)')

# Set limits of plots:
for name in ['visit1', 'visit2']:

    ax[name].set_xlim(1.123, 1.67)
    ax[name].set_ylim(4300,5500)

folder = '../figure3/retrieval-results/hst-retrievals-CLR-FINAL-R10000-POSEIDONBINNING'

# All right --- now we are ready to plot things. Let's do a for loop to generate the plots:
espinoza_color = 'black'#'#6883ba'#'cornflowerblue'
atmosphere_color = 'grey'
ms = 8
ew = 2
alpha = 1.0

nmodels = np.load(folder+'/data_multi_visit_hst_gp_H2-H2atmosphere_george_binned_posterior_contamination1.npy').shape[0]
for name, date in zip(['1', '2'], ['Visit 1', 'Visit 2']):

    # Print date:
    ax['visit'+name].text(1.58,4400,date, fontsize = 30)

    # Read data:
    w_ne, d_ne, derr_ne = np.loadtxt('../figure3/data/hst/visit'+name+'.txt', unpack = True, usecols = (0, 3, 4)) 

    # Extract posterior model samples from both the spectrum and the contamination component:
    contamination = np.load(folder+'/data_multi_visit_hst_gp_H2-H2atmosphere_george_binned_posterior_contamination'+name+'.npy')
    model = np.load(folder+'/data_multi_visit_hst_gp_H2-H2atmosphere_george_binned_posterior_spectrum'+name+'.npy')
    full_model = contamination * model

    # Generate median contamination model for each case:
    median_contamination = np.nanmedian( contamination.T * np.nanmedian( model, axis = 1 ), axis = 1 )
    #median_contamination_fl = np.nanmedian( contamination_fl.T * np.nanmedian( model_fl, axis = 1 ), axis = 1 )

    # Now, get 1, 3 and 5 sigma bands for each full model:
    full_model_median = np.zeros(full_model.shape[1])
    full_model_1sigma_up = np.zeros(full_model.shape[1])
    full_model_1sigma_down = np.zeros(full_model.shape[1])
    full_model_3sigma_up = np.zeros(full_model.shape[1])
    full_model_3sigma_down = np.zeros(full_model.shape[1])

    for i in range(full_model.shape[1]):

        full_model_median[i], full_model_1sigma_up[i], full_model_1sigma_down[i] = juliet.utils.get_quantiles(full_model[:, i])

        _, full_model_3sigma_up[i], full_model_3sigma_down[i] = juliet.utils.get_quantiles(full_model[:, i], alpha = 0.99)

    # Now plot median models for both:
    ax['visit'+name].plot(w_ne, full_model_median, lw = 6, color = atmosphere_color, zorder = 4)
    #ax['visit'+name].plot(w_ne, median_contamination - 250, '--', lw = 3, color = atmosphere_color, zorder = 4)

    # Plot 1 and 3-sigma bands:
    ax['visit'+name].fill_between(w_ne, full_model_1sigma_down, full_model_1sigma_up, color = atmosphere_color, alpha = 0.3, zorder = 2)
    ax['visit'+name].fill_between(w_ne, full_model_3sigma_down, full_model_3sigma_up, color = atmosphere_color, alpha = 0.2, zorder = 2)

    # Plot data:
    ax['visit'+name].errorbar(w_ne, d_ne, derr_ne, fmt = 'o', \
                      ecolor = espinoza_color, \
                      mfc = espinoza_color, \
                      mec = espinoza_color, \
                      ms = ms, \
                      elinewidth = ew, \
                      label = 'NE', alpha = alpha, zorder = 10) 

plt.tight_layout()
plt.savefig('tspectra_bands_'+folder.split('/')[-1]+'.pdf',
            dpi=250,#, rasterize=True,
            #transparent=True,
            bbox_inches='tight')
