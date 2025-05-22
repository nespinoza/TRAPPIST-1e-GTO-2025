import batman, catwoman
import os, sys
import pickle
from datetime import datetime
import numpy as np
import matplotlib
from matplotlib import cm
from matplotlib.patches import Arrow
from matplotlib.gridspec import GridSpec
from scipy.ndimage import gaussian_filter
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

fig = plt.figure(figsize=(15,8))
fig.set_facecolor('w')

# First, cut the parts where the spectra will go:
plt.ylabel('Relative transit depth (ppm)')
plt.xlabel(r'Wavelength ($\mu$m)')
plt.xlim(1.123, 1.67)
plt.ylim(-200,200)

yticks = np.arange(-200,250,50)
plt.yticks(yticks, yticks.astype('str'))

# All right --- now we are ready to plot things. Let's do a for loop to generate the plots:
espinoza_color = 'black'#'#6883ba'#'cornflowerblue'
atmosphere_color = 'grey'#'#1E90FF'
band_atmosphere_line_color = 'grey'
bf_model = 'noflat'

folder = '../figure3/retrieval-results/hst-retrievals-CLR-FINAL-R10000-POSEIDONBINNING'
plot_earth_like = False

if plot_earth_like:

    earth_like_color = 'forestgreen'
    w_el, d_el = np.loadtxt('../figure3/TRAPPIST1e_noerror_earthlike.csv', unpack = True, usecols = (0,1), delimiter = ',')

ms = 13
ew = 4
alpha = 0.1

if 'ghostmu' in folder:

    prefix = 'data_multi_visit_espinoza_gp_ghost-ghostatmosphere'

else:

    prefix = 'data_multi_visit_hst_gp_H2-H2atmosphere'

results = pickle.load(open(folder+'/'+prefix+'_george_results.pkl', 'rb'))
#T, eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO, Pcloud, offset1, rho1, gp_sigma1, sigma_w1, \
#                    offset2, rho2, gp_sigma2, sigma_w2, \
#                    offset3, rho3, gp_sigma3, sigma_w3, \
#                    offset4, rho4, gp_sigma4, sigma_w4 = results['posterior_samples'].T

nmodels = np.load(folder+'/'+prefix+'_george_binned_posterior_contamination1.npy').shape[0]

combined_data = np.array([])
combined_data_err = np.array([])

combined_models = np.array([])
combined_models_fl = np.array([])

counter = 0

for name, date in zip(['1', '2'], ['Visit 1', 'Visit 2']):

    # Read data:
    w_ne, d_ne, derr_ne = np.loadtxt('../figure3/data/hst/visit'+name+'.txt', unpack = True, usecols = (0, 3, 4)) 

    # Extract posterior model samples from both the spectrum and the contamination component:
    contamination = np.load(folder+'/'+prefix+'_george_binned_posterior_contamination'+name+'.npy')
    model = np.load(folder+'/'+prefix+'_george_binned_posterior_spectrum'+name+'.npy')
    full_model = contamination * model

    # Generate median contamination model for each case:
    if bf_model == 'flat':

        median_contamination = np.nanmedian( contamination_fl.T, axis = 1 )

    else:

        median_contamination = np.nanmedian( contamination.T, axis = 1 )

    # Plot decontaminated data:
    corrected_spectra = d_ne / median_contamination

    # Median removal:
    median_value = np.nanmedian( model )#corrected_spectra )

    if counter == 0:

        combined_models = (full_model / median_contamination ) - median_value
        combined_data = corrected_spectra - median_value
        combined_data_err = derr_ne / median_contamination
        counter += 1

    else:

        combined_models = np.vstack(( combined_models, ( full_model / median_contamination ) - median_value))
        combined_data = np.vstack(( combined_data, corrected_spectra - median_value ))
        combined_data_err = np.vstack(( combined_data_err, derr_ne / median_contamination ))

#    plt.errorbar(w_ne, corrected_spectra - np.nanmedian(corrected_spectra), derr_ne / median_contamination, fmt = 'o', \
#                       ecolor = espinoza_color, \
#                       mfc = espinoza_color, \
#                       mec = espinoza_color, \
#                       ms = ms, \
#                       elinewidth = ew, \
#                       label = 'NE', alpha = alpha, zorder = 10) 

# Monte-carlo average data to get errorbars:
combined = np.zeros(combined_data.shape[1])
combined_err = np.zeros(combined_data.shape[1])
for i in range(combined_data.shape[1]):

    avg_distribution = 0. 
    for j in range(combined_data.shape[0]):

        avg_distribution += np.random.normal(combined_data[j,i], combined_data_err[j, i], 1000) 

    avg_distribution = avg_distribution / combined_data.shape[0]

    combined[i] = np.mean(avg_distribution)
    combined_err[i] = np.sqrt( np.var(avg_distribution) )

# Generate models with errorbands:
m = np.zeros(len(w_ne))
m_1sigma_up = np.zeros(len(w_ne))
m_3sigma_up = np.zeros(len(w_ne))
m_1sigma_down = np.zeros(len(w_ne))
m_3sigma_down = np.zeros(len(w_ne))

for i in range(len(w_ne)):

    m[i], m_1sigma_up[i], m_1sigma_down[i] = juliet.utils.get_quantiles(combined_models[:, i])
    _, m_3sigma_up[i], m_3sigma_down[i] = juliet.utils.get_quantiles(combined_models[:, i], alpha = 0.99) 

# Now plot median models for both:
plt.plot(w_ne, m, lw = 7, color = atmosphere_color, zorder = 6)

# Plot EL model if wanted:
if plot_earth_like:

    idx = np.where(w_el < 5.0)[0]
    d_el = gaussian_filter(d_el, 5)
    plt.plot(w_el, d_el - np.nanmedian(d_el[idx]), lw = 7, color = 'grey', zorder = 8)

alphas = 0.2

# Plot 1 and 3-sigma bands:
plt.fill_between(w_ne, m_1sigma_down, m_1sigma_up, color = band_atmosphere_line_color, alpha = alphas, zorder = 2)
plt.fill_between(w_ne, m_3sigma_down, m_3sigma_up, color = band_atmosphere_line_color, alpha = alphas, zorder = 2)

# Plot data:
plt.errorbar(w_ne, combined, combined_err, fmt = 'o',\
                   ecolor = espinoza_color, \
                   mfc = espinoza_color, \
                   mec = espinoza_color, \
                   ms = ms, \
                   elinewidth = ew, \
                   zorder = 11) 

# Print text file:
fout = open('tspectra_decontaminated_corrected_'+bf_model+'_'+folder.split('/')[-1]+'.txt', 'w')
for i in range(len(w_ne)):

    fout.write('{0:.3f} {1:.2f} {2:.2f}\n'.format(w_ne[i], combined[i], combined_err[i]))

fout.close()

plt.tight_layout()
plt.savefig('tspectra_decontaminated_corrected_'+bf_model+'_'+folder.split('/')[-1]+'.pdf',
            dpi=250,#, rasterize=True,
            #transparent=True,
            bbox_inches='tight')
