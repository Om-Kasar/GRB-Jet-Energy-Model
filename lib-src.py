import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Imports the compiled C++ module containing all the functions needed.

# Helper class to load GRB configurations
class BurstConfigurationsLoader:
    def __init__(self, parameters_path, calculated_parameters_path):
        self.parameters = self._load_yaml(parameters_path)
        self.calculated_parameters = self._load_yaml(calculated_parameters_path)
        self.burst_configs = JEF.load_burst_configs({**self.parameters, **self.calculated_parameters})

    @staticmethod
    def _load_yaml(file_path):
        with open(file_path) as file:
            return yaml.safe_load(file)


# Handles jet energy calculations
class JetEnergy:
    def __init__(self, grb_config):
        self.config = grb_config.burst_configs
        self.t_values1 = np.linspace(1.00, self.config.t_NR, 500).tolist()
        self.t_values2 = np.linspace(self.config.t_NR, 1.00e+11, 500).tolist()

    def calculate_jet_energy(self):
        # Calculate the jet energy for two time ranges (before and after t_NR).
        energy1 = JEF.JET_ENERGY1(self.t_values1, self.config)
        energy2 = JEF.JET_ENERGY2(self.t_values2, self.config)
        return energy1, energy2

    def plot(self, ax):
        # Calculate the jet energy
        energy1, energy2 = self.calculate_jet_energy()

        # Plot jet energy
        ax.loglog(self.t_values1, energy1, color='red', label="Jet Energy - Before t_NR")
        ax.loglog(self.t_values2, energy2, color='red', label="Jet Energy - After t_NR")

        # Plot t = t_NR line
        ax.axvline(x=self.config.t_NR, color='red', linestyle='--', label=r"$t = t_{{NR}}$")

        # Formatting the plot
        ax.set_xlim(1.00, 1.00e+11)
        ax.set_ylim(1.00e+43, 1.00e+55)
        ax.set_title('Jet Energy Plot')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Energy [erg]')
        ax.grid(True)
        ax.legend()


# Handles energy per solid angle calculations
class EnergyPerSolidAngleCalculation:
    def __init__(self, grb_config):
        self.config = grb_config.burst_configs
        self.t_values1 = np.linspace(1.00, self.config.t_NR, 1000).tolist()
        self.t_values2 = np.linspace(self.config.t_NR, 1.00e+12, 1000).tolist()

    def calculate_ejection_angles(self):
        """
        Calculate the ejection angles based on time values before and after t_NR.
        """
        theta_rad1 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values1, self.config)]
        theta_rad2 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values2, self.config)]
        return np.degrees(theta_rad1), np.degrees(theta_rad2)

    def calculate_energy_per_solid_angle(self):
        """
        Calculate the energy per solid angle before and after t_NR.
        """
        energies_before_t_NR = JEF.ENERGY_PER_SOLID_ANGLE1(self.t_values1, self.config)
        energies_after_t_NR = JEF.ENERGY_PER_SOLID_ANGLE2(self.t_values2, self.config)
        return energies_before_t_NR, energies_after_t_NR

    def calculate_theta_nr(self):
        """
        Calculate theta_NR in degrees.
        """
        theta_nr_rad = np.array(JEF.ejection_angle([self.config.t_NR], self.config)) - np.array(self.config.initialAngle)
        return np.degrees(theta_nr_rad)

    def plot(self, ax):
        """
        Generates a plot for energy per solid angle with comparison models.
        """
        # Calculate the angle deviations
        change_in_theta_deg1, change_in_theta_deg2 = self.calculate_ejection_angles()

        # Calculate energy per solid angle
        energies_before_t_NR, energies_after_t_NR = self.calculate_energy_per_solid_angle()

        # Calculate theta_NR
        theta_nr_deg = self.calculate_theta_nr()

        # Plot energy per solid angle function
        ax.loglog(
            change_in_theta_deg1,
            energies_before_t_NR,
            color='blue',
            label='Energy per Solid Angle - Our Work'
        )
        ax.loglog(
            change_in_theta_deg2,
            energies_after_t_NR,
            color='blue'
        )

        # Plot theta_NR line
        ax.axvline(x=theta_nr_deg, color='blue', linestyle='--', label=r"$\theta_{NR}$")


# A new class that handles plotting of comparative models (Gaussian and Collapsar)
class ModelPlotter:
    def __init__(self, grb_config):
        self.config = grb_config.burst_configs
        self.t_values1 = np.linspace(1.00, self.config.t_NR, 1000).tolist()
        self.t_values2 = np.linspace(self.config.t_NR, 1.00e+12, 1000).tolist()

    def plot_comparison_models(self, star_mass, ax):
        """
        Plot both the Gaussian and Collapsar models for comparison.
        """
        change_in_theta_deg1, change_in_theta_deg2 = self.calculate_ejection_angles()

        # Plot Gaussian and Collapsar models
        self.plot_gaussian_model(change_in_theta_deg1, color='purple', theta_c=3, ax=ax)
        self.plot_gaussian_model(change_in_theta_deg2, color='orange', theta_c=4, ax=ax)
        self.plot_collapsar_model(change_in_theta_deg1, color='green', star_mass=star_mass, ax=ax)
        self.plot_collapsar_model(change_in_theta_deg2, color='green', star_mass=star_mass, ax=ax)

    def calculate_ejection_angles(self):
        """
        Calculate the ejection angles based on time values before and after t_NR.
        """
        theta_rad1 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values1, self.config)]
        theta_rad2 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values2, self.config)]
        return np.degrees(theta_rad1), np.degrees(theta_rad2)

    def plot_gaussian_model(self, change_in_theta_deg, color, theta_c, ax, label_suffix=""):
        """
        Plots the Gaussian model on the current plot.
        """
        ax.loglog(
            change_in_theta_deg,
            JEF.gaussian_model(change_in_theta_deg, self.config, theta_c),
            color=color,
            label=rf'Gaussian Model ($\theta_c = {theta_c}${label_suffix})'
        )

    def plot_collapsar_model(self, change_in_theta_deg, color, star_mass, ax, label_suffix=""):
        """
        Plots the Collapsar model on the current plot.
        """
        ax.loglog(
            change_in_theta_deg,
            JEF.collapsar_model(change_in_theta_deg, self.config, star_mass),
            color=color,
            label=f'Collapsar Model {label_suffix} (Star Mass = {star_mass} g)'
        )

# Example usage with user-defined star mass
star_mass_input = 1.00e+9  # Change this depending on where the collapsar model (green line) goes in the graph's execution.

# Load GRB2 configurations (modify paths as needed)
burstConfigs = BurstConfigurationsLoader(
    parameters_path="C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml",
    calculated_parameters_path="C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml"
)


# Create a figure for Jet Energy plot
fig1, ax1 = plt.subplots(figsize=(6, 7))

# Create instances of classes for Jet Energy
jet_energy_calculator = JetEnergy(burstConfigs)
jet_energy_calculator.plot(ax1)

# Create a figure for Energy per Solid Angle plot
fig2, ax2 = plt.subplots(figsize=(7, 5))

energy_calculator = EnergyPerSolidAngleCalculation(burstConfigs)
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

model_plotter = ModelPlotter(burstConfigs)
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