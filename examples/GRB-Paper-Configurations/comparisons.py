import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Imports the compiled C++ module containing all the functions needed.

# Helper function to load configurations
def load_configs(parameters_path, calculated_parameters_path):

    with open(parameters_path) as file:
        parameters = yaml.safe_load(file)

    with open(calculated_parameters_path) as file:
        calculated_parameters = yaml.safe_load(file)

    return JEF.load_burst_configs({**parameters, **calculated_parameters})

# Helper function to calculate the angle from jet axis in degrees
def calculate_angle_from_jet_axis(t_values, burst_configs):
    theta_values_rad = JEF.ejection_angle(t_values, burst_configs)
    change_in_theta_rad = [theta_rad - burst_configs.initialAngle for theta_rad in theta_values_rad]
    return np.degrees(change_in_theta_rad)

# Load GRB2 configurations
burstConfigs2 = load_configs(
    "C:/GRB-JET-ENERGY-MODEL/config/GRB2-Configurations/Parameters-GRB2.yaml",
    "C:/GRB-JET-ENERGY-MODEL/config/GRB2-Configurations/CalculatedParameters-GRB2.yaml",
)

# Time values for before and after the non-relativistic phase
t_values1 = np.linspace(1.00, burstConfigs2.t_NR, 500).tolist()
t_values2 = np.linspace(burstConfigs2.t_NR, 1.00e+11, 500).tolist()

# Calculate angle deviations
change_in_theta_deg1 = calculate_angle_from_jet_axis(t_values1, burstConfigs2)
change_in_theta_deg2 = calculate_angle_from_jet_axis(t_values2, burstConfigs2)

# Calculate theta_NR
theta_NR_deg = np.degrees(
    np.array(JEF.ejection_angle([burstConfigs2.t_NR], burstConfigs2)) - burstConfigs2.initialAngle
)

# Initialize the plot
plt.figure(figsize=(7, 5))

# Plot our the
plt.loglog(change_in_theta_deg1, JEF.ENERGY_PER_SOLID_ANGLE1(t_values1, burstConfigs2), color='blue', label='Our Work')
plt.loglog(change_in_theta_deg2, JEF.ENERGY_PER_SOLID_ANGLE2(t_values2, burstConfigs2), color='blue')

# Plot Gaussian Models
plt.loglog(
    change_in_theta_deg1,
    JEF.gaussian_model(change_in_theta_deg1, burstConfigs2, 3),
    color = 'purple',
    label=rf'Gaussian Model ($\theta_c = 3$)'
)
plt.loglog(
    change_in_theta_deg2,
    JEF.gaussian_model(change_in_theta_deg2, burstConfigs2, 4),
    color = 'orange'
)

# Plot Collapsar Model
star_mass = 1.00e+9
plt.loglog(
    change_in_theta_deg1,
    JEF.collapsar_model(change_in_theta_deg1, burstConfigs2, star_mass),
    color='green',
    label='Collapsar Model'
)
plt.loglog(
    change_in_theta_deg2,
    JEF.collapsar_model(change_in_theta_deg2, burstConfigs2, star_mass),
    color='green'
)

# Plot theta_NR line
plt.axvline(x=theta_NR_deg, color='blue', linestyle='--', label=r"$\theta_{NR}$")

# Configure plot limits, labels, and grid
plt.xlim(1.00e-1, 1.00e+4)
plt.ylim(1.00e+40, 1.00e+53)
plt.xlabel(r"$\Delta \theta$ from jet axis [deg]")
plt.ylabel(r"$dE/d\Omega(\theta) \ or \ \epsilon(t)$ [erg/ster]")
plt.grid(True)
plt.legend()
plt.show()