import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF # Imports the compiled C++ module containing all the functions needed.

with open("C:/GRB-JET-ENERGY-MODEL/config/Parameters.yaml") as file: # Replace with your path to the Parameters.yaml file.
    Parameters = yaml.safe_load(file)

with open("C:/GRB-JET-ENERGY-MODEL/config/CalculatedParameters.yaml") as file: # Replace with your path to the CalculatedParameters.yaml file.
    CalculatedParameters = yaml.safe_load(file)

# Defines "burstConfigs" as all of the parameters from the .yaml files which contain all of the burst configuration values.
burstParameters = {**Parameters, **CalculatedParameters}
burstConfigs = JEF.load_burst_configs(burstParameters)

# Define values for time in seconds (s).
t_values1 = np.linspace(1.00, burstConfigs.t_NR, 500).tolist()
t_values2 = np.linspace(burstConfigs.t_NR, 1.00e+10, 500).tolist()

# Values for theta, before the t = t_NR phase.
theta_values_rad1 = JEF.ejection_angle(t_values1, burstConfigs)
change_in_theta_rad1 = [theta_rad - burstConfigs.initialAngle for theta_rad in theta_values_rad1]
change_in_theta_deg1 = np.degrees(change_in_theta_rad1)

# Values for theta, after the t = t_NR phase.
theta_values_rad2 = JEF.ejection_angle(t_values2, burstConfigs)
change_in_theta_rad2 = [theta_rad - burstConfigs.initialAngle for theta_rad in theta_values_rad2]
change_in_theta_deg2 = np.degrees(change_in_theta_rad2)

# Plot the angle in which the graph becomes non-relativistic.
theta_NR_rad = np.array(JEF.ejection_angle([burstConfigs.t_NR], burstConfigs)) - np.array(burstConfigs.initialAngle)

# Plot the values.
plt.figure(figsize=(10, 6))
plt.loglog(change_in_theta_deg1, JEF.gaussian_model(change_in_theta_deg1), color = 'purple')
plt.loglog(change_in_theta_deg1, JEF.ENERGY_PER_SOLID_ANGLE1(t_values1, burstConfigs), color = 'blue')
plt.loglog(change_in_theta_deg2, JEF.ENERGY_PER_SOLID_ANGLE2(t_values2, burstConfigs), color = 'blue')
plt.axvline(x = np.degrees(theta_NR_rad), color = 'red', linestyle = '--', label = r"$\theta_{NR}$")
plt.xlim(1.00e-1, 1.00e+4)
plt.ylim(1.00e+40, 1.00e+53)

plt.xlabel(r"$\Delta \theta$ from jet axis [deg]")
plt.ylabel(r"$dE/d\Omega(\theta)$ [erg/ster]")
plt.grid(True)
plt.legend()
plt.show()