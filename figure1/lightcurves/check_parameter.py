import pickle
import juliet
import numpy as np

fname = 'GPcelerite-white-light-Matern-quadratic-' 
ns = 26948
mean = 0.
mean1 = 0.
mean2 = 0.
pname = 'mflux'

for name, date in zip(['003', '001', '104', '002'], ['June 22', 'June 28', 'July 23', 'October 28']):

    posterior = pickle.load(open(fname+name+'/_dynesty_DNS_posteriors.pkl','rb')) 
    ps = posterior['posterior_samples']

    if pname == 'depth':

        parameter = (ps['p_p1'][:ns]**2)*1e6

    elif pname == 'sigma_w':

        parameter = ps[pname+'_'+name][:ns]

    elif pname == 'GP_sigma' or pname == 'mflux':

        parameter = ps[pname+'_'+name][:ns]*1e6 

    elif pname == 'GP_rho':

        parameter = ps[pname+'_'+name][:ns]*24.

    elif pname == 'ldc':

        q1, q2 = ps['q1_'+name][:ns], ps['q2_'+name][:ns]
        u1, u2 = juliet.utils.reverse_ld_coeffs('quadratic', q1, q2)

        v, vup, vlow = juliet.utils.get_quantiles(u1)
        print(date, ' u1:', v, '+', vup-v, '-', v-vlow)
        mean1 += u1
        mean2 += u2
        v, vup, vlow = juliet.utils.get_quantiles(u2)
        print(date, ' u2:', v, '+', vup-v, '-', v-vlow)

    else:

        parameter = ps[pname][:ns]

    if pname != 'ldc':

        mean = mean + parameter
        v, vup, vlow = juliet.utils.get_quantiles(parameter)
        print(date, ':', v, '+', vup-v, '-', v-vlow)

if pname != 'ldc':

    mean = mean / 4.
    v, vup, vlow = juliet.utils.get_quantiles(mean)
    print('combined: ', v, '+', vup-v, '-', v-vlow)

else:

    mean1 = mean1 / 4.
    mean2 = mean2 / 4.
    v, vup, vlow = juliet.utils.get_quantiles(mean1)
    print('combined, u1: ', v, '+', vup-v, '-', v-vlow)
    v, vup, vlow = juliet.utils.get_quantiles(mean2)
    print('combined, u2: ', v, '+', vup-v, '-', v-vlow)
