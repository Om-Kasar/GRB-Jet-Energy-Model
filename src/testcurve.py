import numpy as np
import matplotlib.pyplot as plt

# Parameters
A_0 = 1.00e51  # Initial energy per solid angle (10^51)
E_final = 1.00e47  # Final energy per solid angle (10^47)
theta_min = 1  # Start of theta range
theta_max = 100  # End of theta range
t_min = 0  # Start time
t_max = 10  # End time
points = 500  # Number of points

# Time-dependent scaling functions
def A(t):
    return A_0 * (1 - t / t_max)  # Linear decay of initial energy with time

def k(t):
    return np.log(A_0 / E_final) / (theta_max - theta_min)  # Decay constant

def theta(t, time):
    return theta_min + (theta_max - theta_min) * (time / t_max)  # Theta as a function of time

# Generate time values
t_values = np.linspace(t_min, t_max, points)

# Generate theta and energy values as functions of time
theta_values = theta(t_values, t_values)
energy_values = [A(t) * np.exp(-k(t) * (theta_value - theta_min)) for t, theta_value in zip(t_values, theta_values)]

# Plot the results
plt.figure(figsize=(8, 6))
plt.plot(theta_values, energy_values, label=r'$E(\theta, t)$', color='blue')
plt.xscale('log')  # Logarithmic theta-axis
plt.yscale('log')  # Logarithmic energy-axis

# Axis labels and limits
plt.xlabel(r'$\theta(t)$ (log scale)')
plt.ylabel(r'Energy per Solid Angle $E(t)$ (log scale)')
plt.title('Exponential Decay of Energy per Solid Angle Over Time')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.ylim(1.00e51, 1.00e47)
plt.legend()
plt.show()
