import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Imports the compiled C++ module containing all the functions needed.

# Helper function to load configurations
def GRB_paths(Parameters_config_path, CalculatedParameters_config_path):

    with open(Parameters_config_path) as file:
        parameters = yaml.safe_load(file)

    with open(CalculatedParameters_config_path) as file:
        calculated_parameters = yaml.safe_load(file)

    return JEF.load_burst_configs({**parameters, **calculated_parameters})

# Helper function to calculate the angle from jet axis in degrees
def calculate_angle_from_jet_axis(t_values, burst_configs):
    theta_values_rad = JEF.ejection_angle(t_values, burst_configs)
    change_in_theta_rad = [theta_rad - burst_configs.initialAngle for theta_rad in theta_values_rad]
    return np.degrees(change_in_theta_rad)

# Load GRB2 configurations
burstConfigs = GRB_paths(
    Parameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml", # Replace with your path to the Parameters.yaml file.
    CalculatedParameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml", # Replace with your path to the CalculatedParameters.yaml file.
)

# Time values for before and after the non-relativistic phase
t_values_before_t_NR = np.linspace(1.00, burstConfigs.t_NR, 500).tolist()
t_values_after_t_NR = np.linspace(burstConfigs.t_NR, 1.00e+11, 500).tolist()

# Calculate angle deviations
change_in_theta_deg1 = calculate_angle_from_jet_axis(t_values_before_t_NR, burstConfigs)
change_in_theta_deg2 = calculate_angle_from_jet_axis(t_values_after_t_NR, burstConfigs)

# Calculate theta_NR
theta_NR_deg = np.degrees(
    np.array(JEF.ejection_angle([burstConfigs.t_NR], burstConfigs)) - burstConfigs.initialAngle
)

# Initialize the plot
plt.figure(figsize=(7, 5))

# Plot our energy per solid angle function.
plt.loglog(change_in_theta_deg1, JEF.ENERGY_PER_SOLID_ANGLE1(t_values_before_t_NR, burstConfigs), color='blue', label='Our Work')
plt.loglog(change_in_theta_deg2, JEF.ENERGY_PER_SOLID_ANGLE2(t_values_after_t_NR, burstConfigs), color='blue')

# Plot Gaussian Models
plt.loglog(
    change_in_theta_deg1,
    JEF.gaussian_model(change_in_theta_deg1, burstConfigs, 3),
    color = 'purple',
    label=rf'Gaussian Model ($\theta_c = 3$)'
)
plt.loglog(
    change_in_theta_deg2,
    JEF.gaussian_model(change_in_theta_deg2, burstConfigs, 4),
    color = 'orange'
)

# Plot Collapsar Model
star_mass = 1.00e+9
plt.loglog(
    change_in_theta_deg1,
    JEF.collapsar_model(change_in_theta_deg1, burstConfigs, star_mass),
    color='green',
    label='Collapsar Model'
)
plt.loglog(
    change_in_theta_deg2,
    JEF.collapsar_model(change_in_theta_deg2, burstConfigs, star_mass),
    color='green'
)

# Plot theta_NR line
plt.axvline(x = theta_NR_deg, color='blue', linestyle='--', label=r"$\theta_{NR}$")

# Configure plot limits, labels, and grid
plt.xlim(1.00e-1, 1.00e+4)
plt.ylim(1.00e+40, 1.00e+53)
plt.xlabel(r"$\Delta \theta$ from jet axis [deg]")
plt.ylabel(r"$dE/d\Omega(\theta) \ or \ \epsilon(t)$ [erg/ster]")
plt.grid(True)
plt.legend()
plt.show()