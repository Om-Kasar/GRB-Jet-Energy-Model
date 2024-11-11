import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import computationalfunctions as cf

t_values = np.linspace(0.0001, 1e+9, 1000)

lorentz_factor_values = [cf.Lorentz_Factor(t_val) for t_val in t_values]

plt.plot(t_values, lorentz_factor_values, color = 'green')
plt.show()