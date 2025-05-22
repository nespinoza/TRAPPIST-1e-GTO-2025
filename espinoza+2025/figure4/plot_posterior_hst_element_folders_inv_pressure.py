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
import arviz as az

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


def eps_to_MR(eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO):
 
         # First, compute the eps of the reminder by the condition sum(eps) 
         eps_H2 = - (eps_CO2 + eps_CH4 + eps_H2O + eps_N2 + eps_O2 + eps_O3 + eps_N2O + eps_CO)
     
         # Next, from all the eps values, compute the value of the g(X) funct
         g = 1. / (np.exp(eps_H2) + np.exp(eps_CO2) + np.exp(eps_CH4) + np.exp(eps_H2O) + np.exp(eps_N2) + np.exp(eps_O2)
     + np.exp(eps_O3) + np.exp(eps_N2O) + np.exp(eps_CO))
     
         # With this, get back the mixing ratios for the eleme
         return np.exp(eps_H2) * g, np.exp(eps_CO2) * g, np.exp(eps_CH4) * g, np.exp(eps_H2O) * g, np.exp(eps_N2) * g, np.exp(eps_O2) * g, np.exp(eps_O3) * g, np.exp(eps_N2O) * g, np.exp(eps_CO) * g
    

# Get JWST outp
element = 'H2'
folder = 'hst-retrievals-CLR-FINAL-R10000-POSEIDONBINNING'
results_JWST = pickle.load(open('../figure3/retrieval-results/'+folder+'/data_multi_visit_hst_gp_H2-H2atmosphere_george_results.pkl', 'rb'))

T, eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO, jwst_Pcloud, offset1, rho1, gp_sigma1, sigma_w1, \
   offset2, rho2, gp_sigma2, sigma_w2, P_ref = results_JWST['posterior_samples'].T

jwst_H2, CO2, CH4, H2O, N2, O2, O3, N2O, CO = eps_to_MR(eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO)
# some random data


y = np.log10(jwst_Pcloud)#np.random.randn(1000)
if element == 'CH4':

    x = np.log10(CH4)#np.random.randn(1000)

elif element == 'H2':

    x = np.log10(jwst_H2)

idx = np.where(y<0.5)[0]
print('posterior mass below 0.5:',(len(idx)/len(y))*100)

pltparams = load_plt_params()

fig = plt.figure(figsize=(8,8))
fig.set_facecolor('w')
modelcolor = 'black'

ax = fig.add_subplot(111)

az.plot_kde(
x,
y,
hdi_probs=[0.393, 0.865, 0.989],  # 1, 2 and 3 sigma contours
contourf_kwargs={"cmap": "Blues"},
ax = ax
)

# value for Mars from https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021EA001915 (max value)
#ax.plot([np.log10(30*1e-9), np.log10(30*1e-9)], [-6, 2.], '--', color = 'red')
# Same for venus from https://www.sciencedirect.com/science/article/pii/S2215016123003023
#ax.plot([np.log10(300*1e-6), np.log10(300*1e-6)], [-6, 2.], '--', color = 'orange')

#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel(r'$\log$ P$_{surf}$ (bar)', fontsize = 25)
ax.set_xlabel('$\log$ '+element, fontsize = 25)
ax.set_ylim(2., -6)
ax.set_xlim(-12.,0.)

ax.set_yticks([-6,-4,-2, 0], ['-6', '-4', '-2', '0'])
ax.set_xticks([-10, -7.5, -5.0, -2.5], ['$-10$', '$-7.5$', '$-5.0$', '$-2.5$'])

print('saving to','posteriors_hst_log_'+element+'_'+folder+'.pdf')
plt.savefig('posteriors_hst_log_'+element+'_'+folder+'-ip.pdf',
            dpi=250,#, rasterize=True,
            #transparent=True,
            bbox_inches='tight')
