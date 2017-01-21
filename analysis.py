import data as d
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

Cu0 = d.Data("Ni_11_1_17.txt").lattice_parameter()
Cu25 = d.Data("Cu_25_Ni_75_11_1_17.txt").lattice_parameter()
Cu50 = d.Data("Cu50Ni50_11_1_17.txt").lattice_parameter()
Cu75 = d.Data("Cu_75_Ni_25_12_1_17.txt").lattice_parameter()
Cu100 = d.Data("Cu_11-01-17.txt").lattice_parameter()

x = np.array([0, 25, 50, 75, 100])
y = np.array([Cu0[0], Cu25[0], Cu50[0], Cu75[0], Cu100[0]])
err = np.array([Cu0[1], Cu25[1], Cu50[1], Cu75[1], Cu100[1]])

def linear(x, a, b):
    return a * x + b

popt, pcov = curve_fit(
            linear,
            x, y,
            sigma=err,
            p0=(0.1/100.0, 3.6)
        )

y_fit = popt[0]*x + popt[1]
chi_square = (y - y_fit)**2 / err**2
print chi_square

plt.errorbar(x, y, yerr=err)
plt.plot(x, y_fit)
plt.show()
