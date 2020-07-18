# Script to run Challenge 19-2

from infection_model import infection_model

m = 10
n = 10
k = 4
tau = 0.2
delta = 0.01

model = infection_model(m = m, n = n, k = k, tau = tau, delta = delta)
print(model.run_simulation(True))
model.plot_data(save = True, show = True)