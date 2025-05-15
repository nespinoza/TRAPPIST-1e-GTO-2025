import juliet
import matplotlib.pyplot as plt
import numpy as np
import pickle

ns=26948 
posteriors22 = pickle.load(open('GPcelerite-white-light-Matern-quadratic-003/_dynesty_DNS_posteriors.pkl','rb'))
posteriors28 = pickle.load(open('GPcelerite-white-light-Matern-quadratic-001/_dynesty_DNS_posteriors.pkl','rb'))

d28 = (posteriors28['posterior_samples']['p_p1']**2)*1e6
d22 = (posteriors22['posterior_samples']['p_p1']**2)*1e6

d22 = d22[:ns]
d28 = d28[:ns]

plt.hist(d22, label ='22 june')
plt.hist(d28, label = '28 june')

plt.show()

plt.hist(d22-d28, bins = 100)
difference = d22-d28
v,vup,vlow = juliet.utils.get_quantiles(difference)
print(v, '+',vup-v, '-',v-vlow)
plt.show()
