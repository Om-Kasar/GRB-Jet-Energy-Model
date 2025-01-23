import jet_energy_functions as jef
import matplotlib.pyplot as plt

star_mass_input = 1.00e+9  # Collapsar Model

# Load GRB2 configurations (modify paths as needed)
burstConfigs1 = jef.load_burst_configs(
    parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB1-Configurations\\Parameters-GRB1.yaml",
    calculated_parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB1-Configurations\\CalculatedParameters-GRB1.yaml"
)

burstConfigs2 = jef.load_burst_configs(
    parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB2-Configurations\\Parameters-GRB2.yaml",
    calculated_parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB2-Configurations\\CalculatedParameters-GRB2.yaml"
)

burstConfigs3 = jef.load_burst_configs(
    parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB3-Configurations\\Parameters-GRB3.yaml",
    calculated_parameters_path="C:\\repos\\GRB-Jet-Energy-Model\\examples\\GRB-Paper-Configurations\\config\\GRB3-Configurations\\CalculatedParameters-GRB3.yaml"
)

# Create a figure for Jet Energy plot
fig1, ax1 = plt.subplots(figsize=(6, 7))

# Create instances of classes for Jet Energy
jet_energy_calculator1 = jef.jet_energy(burstConfigs1)
jet_energy_calculator2 = jef.jet_energy(burstConfigs2)
jet_energy_calculator3 = jef.jet_energy(burstConfigs3)

jet_energy_calculator1.plot(ax1, color='red', label=1)
jet_energy_calculator2.plot(ax1, color='blue', label=2)
jet_energy_calculator3.plot(ax1, color='orange', label=3)

# Create a figure for Energy per Solid Angle plot
fig2, ax2 = plt.subplots(figsize=(7, 5))

energy_calculator1 = jef.energy_per_solid_angle(burstConfigs1)
energy_calculator2 = jef.energy_per_solid_angle(burstConfigs2)
energy_calculator3 = jef.energy_per_solid_angle(burstConfigs3)

energy_calculator1.plot(ax2, color='red', label=1)
energy_calculator2.plot(ax2, color='blue', label=2)
energy_calculator3.plot(ax2, color='orange', label=3)

# Formatting the plot for Energy per Solid Angle
ax2.set_xlim(1.00e-1, 1.00e+4)
ax2.set_ylim(1.00e+40, 1.00e+53)
ax2.set_title('Energy per Solid Angle')
ax2.set_xlabel(r"$\Delta \theta$ from jet axis [deg]")
ax2.set_ylabel(r"$dE/d\Omega(\theta)$ or $\epsilon(t)$ [erg/ster]")
ax2.grid(True)
ax2.legend()

# Create a figure for Energy per Solid Angle and Comparison Models plot
fig3, ax3 = plt.subplots(figsize=(7, 5))

model_plotter = jef.model_plotter(burstConfigs2)
energy_calculator2.plot(ax3, color='blue', label=2)
model_plotter.plot_comparison_models(star_mass=star_mass_input, ax=ax3)

# Formatting the plot for Energy per Solid Angle and Comparison Models
ax3.set_xlim(1.00e-1, 1.00e+4)
ax3.set_ylim(1.00e+40, 1.00e+54)
ax3.set_title('Energy per Solid Angle and Comparative Models')
ax3.set_xlabel(r"$\Delta \theta$ from jet axis [deg]")
ax3.set_ylabel(r"$dE/d\Omega(\theta)$ or $\epsilon(t)$ [erg/ster]")
ax3.grid(True)
ax3.legend()

# Show all plots at once
plt.show()