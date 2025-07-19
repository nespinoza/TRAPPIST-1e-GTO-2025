import numpy as np
from scipy.stats import chi2

for i in [1, 2, 3, 4]:

    print('visit',i)
    print('-----------')
    w,d,derr = np.loadtxt('visit'+str(i)+'.txt', unpack = True, usecols = (0, 3, 4))
    ws,ds,derrs = np.loadtxt('stevenson/visit'+str(i)+'.txt', unpack = True, usecols = (0, 3, 4))
    wg,dg,derrg = np.loadtxt('grant/visit'+str(i)+'.txt', unpack = True, usecols = (0, 3, 4))

    d = d-np.nanmedian(d)
    ds = ds-np.nanmedian(ds)
    dg = dg-np.nanmedian(dg)

    residuals = d - ds
    residual_err = np.sqrt(derr**2 + derrs**2)
    chi = residuals / residual_err

    chi2_stat = np.sum(chi**2)
    dof = len(residuals)
    p_value = chi2.sf(chi2_stat, dof)

    print('stev:')
    print(p_value)

    print('grant:')
    residuals = d[3:] - dg
    residual_err = np.sqrt(derr[3:]**2 + derrg**2)
    chi = residuals / residual_err
    chi2_stat = np.sum(chi**2)
    dof = len(residuals)
    p_value = chi2.sf(chi2_stat, dof)
    print(p_value)
    print('\n')
