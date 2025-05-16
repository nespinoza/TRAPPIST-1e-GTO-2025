import batman, catwoman
import os, sys
import pickle
from datetime import datetime
import numpy as np
import matplotlib
from scipy.ndimage import median_filter
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

fig = plt.figure(figsize=(25,9))
#fig = plt.figure(figsize=(25,6))
#fig = plt.figure(figsize=(26,6))
fig.set_facecolor('w')
modelcolor = 'black'

# First, cut the parts where the lightcurves will go:
height = 330
gs = GridSpec(height, 4, figure=fig, wspace = 0.1)
ax1 = {}
ax2 = {}
ax3 = {}

transit_height = 140
delta_plot = 20
no_transit_height = height - transit_height - delta_plot * 2

print(transit_height+int(no_transit_height*0.5)+delta_plot*2)
for name, index in zip(['003', '001', '104', '002'], [0, 1, 2, 3]):

    ax1[name] = fig.add_subplot(gs[:transit_height, index])
    ax2[name] = fig.add_subplot(gs[transit_height+delta_plot:transit_height+delta_plot+int(no_transit_height*0.5), index])
    ax3[name] = fig.add_subplot(gs[transit_height+int(no_transit_height*0.5)+delta_plot*2:, index])

# Set labels:
ax1['003'].set_ylabel('Relative flux')
ax2['003'].set_ylabel('O-C (ppm)')
ax3['003'].set_ylabel(r'H$\alpha$ (ppt)')

for name in ['003', '001', '104', '002']:

    # Set h-alpha plot x-labels:
    ax3[name].set_xlabel('Time from mid-transit (hr)')
    
    # Take the chance to remove ticks from top plots in the bottom:
    ax1[name].tick_params(labelbottom=False)
    ax2[name].tick_params(labelbottom=False)

# Remove ticks from y-axes of all plots to the right of the first:
for name in ['001', '104', '002']:

    ax1[name].tick_params(labelleft=False)
    ax2[name].tick_params(labelleft=False)
    ax3[name].tick_params(labelleft=False)

# Set limits of plots:
for name in ['003', '001', '104', '002']:

    ax1[name].set_xlim(-0.8, 0.8)
    ax2[name].set_xlim(-0.8, 0.8)
    ax3[name].set_xlim(-0.8, 0.8)

    ax1[name].set_ylim(0.993, 1.0015)
    ax2[name].set_ylim(-950,950)
    ax3[name].set_ylim(-70,70)

# All right --- now we are ready to plot things. Let's do a for loop to generate the plots:
for name, date in zip(['003', '001', '104', '002'], ['June 22', 'June 28', 'July 23', 'October 28']):

    dataset = juliet.load(input_folder = 'lightcurves/GPcelerite-white-light-Matern-quadratic-'+name)
    results = dataset.fit(sampler = 'dynamic_dynesty')
    model = results.lc.evaluate(name)

    # Extract time-of-transit center and jitter term:
    t0 = np.nanmedian( results.posteriors['posterior_samples']['t0_p1'] )
    sigma_w = np.nanmedian( results.posteriors['posterior_samples']['sigma_w_'+name] )

    # Extract times and fluxes:
    t, f, ferr = dataset.times_lc[name], dataset.data_lc[name], dataset.errors_lc[name]
    # Add jitter term to errorbars:
    ferr = np.sqrt( ferr**2 + (sigma_w*1e-6)**2 )

    thour = ( t - t0 ) * 24

    # Generate top plot:
    ax1[name].set_title(date+', 2023', fontsize = 25)

    ax1[name].errorbar( thour, f, ferr, fmt = '.', ms = 1, elinewidth = 1, color = 'black', alpha = 0.3 )
    ax1[name].plot( thour, model, color = modelcolor, lw = 3 )

    # Get bottom/residual plot:
    r = (f - model) * 1e6
    ax2[name].errorbar( thour, r, ferr*1e6, fmt = '.', ms = 1, elinewidth = 1, color = 'black', alpha = 0.3 )
    ax2[name].plot( [-10,10], [0, 0], '--', color = modelcolor, lw = 3 )
    tb = ax2[name].text( 0.0, 550, r'$\sigma = {0:}$ ppm'.format(str(int(np.sqrt(np.var(r))))), fontsize = 25)

    # Extract h-alpha lightcurves:
    month, day = date.split()
    tha, fha, ferrha = np.loadtxt('halpha-lcs/'+month.lower()+'-'+day+'_halpha.txt', unpack = True, usecols = (0,1,2))
    ax3[name].errorbar( tha, ( fha - 1. ) * 1e3, ferrha * 1e3, fmt = '.', color = 'cornflowerblue', ms = 5, elinewidth = 1, alpha = 0.2 )
    tbin, fbin, fbinerr = juliet.utils.bin_data( tha, ( fha - 1. ) * 1e3, yerr = ferrha * 1e3, n_bin = 20 )
    ax3[name].errorbar( tbin, fbin, fbinerr, fmt = 'o', color = 'cornflowerblue', ms = 6, elinewidth = 2)
    #ax3[name].plot( tha, median_filter( ( fha - 1. ) * 1e3, 10 ), lw = 3 )

plt.savefig('lcs.pdf',
            dpi=250,#, rasterize=True,
            #transparent=True,
            bbox_inches='tight')
