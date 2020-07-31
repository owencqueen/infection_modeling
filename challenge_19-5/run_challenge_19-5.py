# Script to run Challenge 19-5
#   Note: Make sure to include infection_model.py file in this directory

import numpy as np
from matplotlib import pyplot as plt
from polyFit import *
from infection_model import infection_model

m = 10
n = 10
k = 4
tau = 0.2
delta = 0.01
nu = 0.1

model = infection_model(m = m, n = n, k = k, tau = tau, nu = nu)

'''
# CODE BLOCK #1
#..........................................................................................
# In order to run the process for finding the nu value at 20% containment, 
#   uncomment this block of code (remove the ''' ''' surrounding this block).
#..........................................................................................

# Run our model for 20 runs with given parameters, take average best nu
# Run with AIC optimization
new_nu = 0
for i in range(0, 20):
    new_nu += model.optimize_nu(min_fn = AIC, runs = 100, nu_0 = 0.0, nu_f = 1, increment = 0.01,
                             plot = False, pct_containment = 0.2)
new_nu /= 20
print("Optimized nu value with 0.2 containment = {:.4f}".format(new_nu))
'''

'''
#  CODE BLOCK #2
#..........................................................................................
# In order to run the process for determining the nu values at each percent of containment, 
#   uncomment this block of code (remove the ''' ''' surrounding this block).
#..........................................................................................

best_nu = []
pct_containment_vals = np.arange(start = 0.02, stop = 0.21, step = 0.01)
for pct in pct_containment_vals:
    bnu = model.optimize_nu(min_fn = AIC, runs = 100, nu_0 = 0.0, nu_f = 1, increment = 0.01,
                             plot = False, pct_containment = pct)
    print('Calculated nu value with {:.2f} containment = {:.4f}'.format(pct, bnu))
    best_nu.append(bnu)

plt.plot(pct_containment_vals, best_nu, "-", color = "red")

plt.title("Nu Needed for given Containment Percentage")
plt.xlabel("Percent Containment")
plt.ylabel("Nu")

plt.show()
'''