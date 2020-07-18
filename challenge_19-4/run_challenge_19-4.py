# Script to run Challenge 19-4

from infection_model import infection_model

m = 10
n = 10
k = 4
tau = 0.2
delta = 0.01
nu = 0.1

model = infection_model(m = m, n = n, k = k, tau = tau, nu = nu)
model.run_multiple_sim(runs = 100, nu_0 = 0.0, nu_f = 0.3, increment = 0.1)
model.plot_histograms_infected(save = False, show = True)
