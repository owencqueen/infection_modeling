# Script to run Challenge 19-3
#   Note: Make sure to include infection_model.py file in this directory

from infection_model import infection_model

m = 10
n = 10
k = 4
tau = 0.2
delta = 0.01
nu = 0.1

model = infection_model(m = m, n = n, k = k, tau = tau, nu = nu)
print(model.run_simulation(True))
model.plot_data(plot_title = "Infection Model with Vaccination", save = True, show = True)