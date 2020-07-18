# Script to run challenge 19-1

from infection_model import *

m = 10
n = 10
k = 4
tau = 0.2

model = infection_model(m = m, n = n, k = k, tau = tau)
print(model.run_simulation(True))
model.plot_data(save = True, show = True)