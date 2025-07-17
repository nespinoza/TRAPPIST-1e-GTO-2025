import numpy as np
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

import corner

def eps_to_MR(eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO):
 
         # First, compute the eps of the reminder by the condition sum(eps) 
         eps_H2 = - (eps_CO2 + eps_CH4 + eps_H2O + eps_N2 + eps_O2 + eps_O3 + eps_N2O + eps_CO)
    
         # Next, from all the eps values, compute the value of the g(X) funct
         g = 1. / (np.exp(eps_H2) + np.exp(eps_CO2) + np.exp(eps_CH4) + np.exp(eps_H2O) + np.exp(eps_N2) + np.exp(eps_O2)
     + np.exp(eps_O3) + np.exp(eps_N2O) + np.exp(eps_CO))
    
         # With this, get back the mixing ratios for the eleme
         return np.exp(eps_H2) * g, np.exp(eps_CO2) * g, np.exp(eps_CH4) * g, np.exp(eps_H2O) * g, np.exp(eps_N2) * g, np.exp(eps_O2) * g, np.exp(eps_O3) * g, np.exp(eps_N2O) * g, np.exp(eps_CO) * g


results_JWST = pickle.load(open('../figure3/retrieval-results/JWST-CLR-R10000/data_multi_visit_espinoza_gp_H2-H2atmosphere_george_results.pkl','rb'))

T, eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO, jwst_Pcloud, offset1, rho1, gp_sigma1, sigma_w1, \
   offset2, rho2, gp_sigma2, sigma_w2, \
   offset3, rho3, gp_sigma3, sigma_w3, \
   offset4, rho4, gp_sigma4, sigma_w4, P_ref = results_JWST['posterior_samples'].T

jwst_H2, CO2, CH4, H2O, N2, O2, O3, N2O, CO = eps_to_MR(eps_CO2, eps_CH4, eps_H2O, eps_N2, eps_O2, eps_O3, eps_N2O, eps_CO)

# To avoid too much data crowding, sample only 10,000 points:
idx = np.random.choice(np.arange(len(jwst_H2)), 10000, replace = False)

data = np.vstack(( np.log10(jwst_H2[idx]), np.log10(jwst_Pcloud[idx]) ))

for param in [rho1, rho2, rho3, rho4, gp_sigma1, gp_sigma2, gp_sigma3, gp_sigma4]:

    data = np.vstack((data, param[idx]))

data = data.T

labels =[ r"$\log H_2$",
          r"$\log P_{cloud}$",
          r"$\ell_1$",
          r"$\ell_2$",
          r"$\ell_3$",
          r"$\ell_4$",
          r"$A_{1}$",
          r"$A_{2}$",
          r"$A_{3}$",
          r"$A_{4}$"
       ]

ranges = [
          (np.min(np.log10(jwst_H2)), np.max(np.log10(jwst_H2))),
          (-6, 2),
          (0,100),
          (0,100),
          (0,100),
          (0, 4),
          (0,10),
          (0,10),
          (0,10),
          (0, 0.05)
]

# Plot it.
figure = corner.corner(
    data,
    labels=labels,
    range = ranges,
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 12},
)

plt.savefig('cplot_raw.pdf')
