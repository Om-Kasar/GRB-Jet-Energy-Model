import numpy as np
import matplotlib.pyplot as plt
import sympy as sm

def DegreesToRadians(degreeMeasure):
    return degreeMeasure * (np.pi / 180)

t_values = np.linspace(0.1, 100000, 1000)

t_dec = 10
alpha = 0.35
Gamma_0 = 2000
initialAngle = 47.0101
radiusConstant = 2e5
powerLawConstant = 1e4

x = sm.symbols("x")

# Define the Lorentz factor with a smooth transition using hyperbolic tangent (tanh) function
def lorentz_factor_tanh(t):
    return Gamma_0 * (t_dec / t)**0.5 * (1 + np.tanh(alpha * (t - t_dec))) / 2 + \
           Gamma_0 * (1 - np.tanh(alpha * (t - t_dec))) / 2

def ejection_angle(t):
    return initialAngle + (initialAngle / lorentz_factor_tanh(t))

# Generate Lorentz factor values for plotting
Gamma_tanh_values = lorentz_factor_tanh(t_values)
theta_tanh_values = ejection_angle(t_values)

def radius(t):
    return radiusConstant * (t ** 2/5)

def powerLawModel(t):
    return powerLawConstant * (radius(t) ** 2/5)

def massSweptUpDensity(t):
    return 4 * powerLawModel(t) * lorentz_factor_tanh(t)



