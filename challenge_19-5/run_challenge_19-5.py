from polyFit import *
from infection_model import infection_model

m = 10
n = 10
k = 4
tau = 0.2
delta = 0.01
nu = 0.1

model = infection_model(m = m, n = n, k = k, tau = tau, nu = nu)
new_nu = model.optimize_nu(min_fn = r2ADJ, runs = 100, nu_0 = 0.0, nu_f = 1, increment = 0.01,
                             plot = True, pct_containment = 0.2)
print(new_nu)