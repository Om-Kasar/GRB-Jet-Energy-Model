import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Import the compiled C++ module

class GRB_paths:
    # Handles loading and storing GRB configuration parameters.

    def __init__(self, Parameters_config_path, CalculatedParameters_config_path):
        self.parameters = self._load_yaml(Parameters_config_path)
        self.calculated_parameters = self._load_yaml(CalculatedParameters_config_path)
        self.burst_configs = JEF.load_burst_configs({**self.parameters, **self.calculated_parameters})

    @staticmethod
    def _load_yaml(file_path):
        with open(file_path) as file:
            return yaml.safe_load(file)

class energy_per_solid_angle_calculator:
    # Handles all calculations for a single GRB.

    def __init__(self, grb_config):
        self.config = grb_config.burst_configs
        self.t_values_before_t_NR = np.linspace(1.00, self.config.t_NR, 1000).tolist()
        self.t_values_after_t_NR = np.linspace(self.config.t_NR, 1.00e+12, 1000).tolist()

    def calculate_ejection_angles(self):
        theta_rad1 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values_before_t_NR, self.config)]
        theta_rad2 = [theta - self.config.initialAngle for theta in JEF.ejection_angle(self.t_values_after_t_NR, self.config)]
        return np.degrees(theta_rad1), np.degrees(theta_rad2)

    def calculate_energy_per_solid_angle(self):
        energies_before_t_NR = JEF.ENERGY_PER_SOLID_ANGLE1(self.t_values_before_t_NR, self.config)
        energies_after_t_NR = JEF.ENERGY_PER_SOLID_ANGLE2(self.t_values_after_t_NR, self.config)
        return energies_before_t_NR, energies_after_t_NR

    def calculate_theta_nr(self):
        theta_nr_rad = np.array(JEF.ejection_angle([self.config.t_NR], self.config)) - np.array(self.config.initialAngle)
        return np.degrees(theta_nr_rad)


class GRBPlotter:
    # Handles energy per solid angle plotting for GRBs.

    @staticmethod
    def plot_grbs(grb_calculators, labels, colors):
        plt.figure(figsize=(6, 7))

        for grb, label, color in zip(grb_calculators, labels, colors):
            change_in_theta_deg1, change_in_theta_deg2 = grb.calculate_ejection_angles()
            energies_before_t_NR, energies_after_t_NR = grb.calculate_energy_per_solid_angle()
            theta_nr = grb.calculate_theta_nr()

            # Plot energy per solid angle.
            plt.loglog(change_in_theta_deg1, energies_before_t_NR, color=color, label=rf"{label} - $\epsilon(t)$")
            plt.loglog(change_in_theta_deg2, energies_after_t_NR, color=color)

            # Plot theta = theta_NR.
            plt.axvline(x=theta_nr, color=color, linestyle='--', label=f"$\\theta_{{NR, {label[-1]}}}$")

        # Formatting the plot
        plt.xlim(1.00e-1, 1.00e+4)
        plt.ylim(1.00e+40, 1.00e+54)
        plt.xlabel(r"$\Delta \theta(t)$ from jet axis [deg]")
        plt.ylabel(r"$dE/d\Omega(\theta)$ or $\epsilon(t)$ [erg/ster]")
        plt.grid(True)
        plt.legend()
        plt.show()

# Main script starts here.
# -------------------------------------------------------------------

# GRB Configurations & plotting.
burstConfigs = GRB_paths(
    # Replace with your path to the Parameters.yaml file.
    Parameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml", 

    # Replace with your path to the CalculatedParameters.yaml file.
    CalculatedParameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml"
)

# GRB Calculations
plotted_grb = energy_per_solid_angle_calculator(burstConfigs)

# Plotting the GRB.
GRBPlotter.plot_grbs(
    [plotted_grb],
    labels=["GRB1"],
    colors=["red"],
)