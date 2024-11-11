from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
import matplotlib.pyplot as plt

# Step 1: Load your light curve data
file_path = r"C:\POLYGENCEPROJECT\GRB200806A.lc"
time, flux = np.loadtxt(file_path, unpack=True)

# Step 2: Polynomial fit of the light curve
degree = 10  # Adjust based on data behavior
coeffs = np.polyfit(time, flux, degree)
poly_func = np.poly1d(coeffs)

# Step 3: Define a function to calculate the integral of the polynomial fit
def energy_emitted(t):
    integral, _ = quad(poly_func, 0, t)  # Integrate from 0 to t
    return integral

# Step 4: Calculate the energy over a range of times
t_values = np.linspace(0, max(time), 100)  # Define multiple t values
energies = [energy_emitted(t) for t in t_values]

# Output or plot the results
plt.plot(t_values, energies, label="Energy Emitted")
plt.xlabel("Time")
plt.ylabel("Total Energy Emitted")
plt.legend()
plt.show()