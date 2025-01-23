import jet_energy_functions as jef
import matplotlib.pyplot as plt

# Example usage with user-defined star mass
star_mass_input = 1.00e+9  # Change this depending on where the collapsar model (green line) goes in the graph's execution.

# Load GRB configurations (modify paths as needed)
burstConfigs = jef.load_burst_configs(
    parameters_path="C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml", # Path to the "Parameters.yaml" file located in this directory.
    calculated_parameters_path="C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml" # Path to the "CalculatedParameters.yaml" file located in this directory.
)

# Create a figure for Jet Energy plot
fig1, ax1 = plt.subplots(figsize=(6, 7))

# Create instances of classes for Jet Energy
jet_energy_calculator = jef.jet_energy(burstConfigs)
jet_energy_calculator.plot(ax1)

# Create a figure for Energy per Solid Angle plot
fig2, ax2 = plt.subplots(figsize=(7, 5))

energy_calculator = jef.energy_per_solid_angle(burstConfigs)
energy_calculator.plot(ax2)

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

model_plotter = jef.model_plotter(burstConfigs)
energy_calculator.plot(ax3)
model_plotter.plot_comparison_models(star_mass=star_mass_input, ax=ax3)

# Formatting the plot for Energy per Solid Angle and Comparison Models
ax3.set_xlim(1.00e-1, 1.00e+4)
ax3.set_ylim(1.00e+40, 1.00e+53)
ax3.set_title('Energy per Solid Angle and Comparative Models')
ax3.set_xlabel(r"$\Delta \theta$ from jet axis [deg]")
ax3.set_ylabel(r"$dE/d\Omega(\theta)$ or $\epsilon(t)$ [erg/ster]")
ax3.grid(True)
ax3.legend()

# Show all plots at once
plt.show()