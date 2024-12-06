import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Imports the compiled C++ module containing all the functions needed for analysis.

class GRB_paths:
    def __init__(self, name, Parameters_config_path, CalculatedParameters_config_path):
        self.name = name
        self.Parameters_config_path = Parameters_config_path
        self.CalculatedParameters_config_path = CalculatedParameters_config_path
        self.parameters = self._load_yaml(f"{Parameters_config_path}")
        self.calculated_parameters = self._load_yaml(f"{CalculatedParameters_config_path}")
        self.burst_configs = JEF.load_burst_configs({**self.parameters, **self.calculated_parameters})
        self.t_values_before_t_NR = np.linspace(1.00, self.burst_configs.t_NR, 500).tolist()
        self.t_values_after_t_NR = np.linspace(self.burst_configs.t_NR, 1.00e+11, 500).tolist()

    def _load_yaml(self, file_path):
        with open(file_path) as file:
            return yaml.safe_load(file)

    def plot(self, ax, color):
        # Plot jet energy for time ranges.
        ax.loglog(self.t_values_before_t_NR, JEF.JET_ENERGY1(self.t_values_before_t_NR, self.burst_configs), color=color, label=f"{self.name}")
        ax.loglog(self.t_values_after_t_NR, JEF.JET_ENERGY2(self.t_values_after_t_NR, self.burst_configs), color=color)
        # Plot t = t_NR line.
        ax.axvline(x=self.burst_configs.t_NR, color=color, linestyle='--', label=fr"$t = t_{{NR, {self.name[-1]}}}$")

# Define the GRB.
burstConfigs = GRB_paths(name="GRB1", 
    # Replace with your path to the Parameters.yaml file.
    Parameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml", 
    # Replace with your path to the CalculatedParameters.yaml file.
    CalculatedParameters_config_path="C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml")

# Plot the GRB.
fig, ax = plt.subplots(figsize=(6, 7))
burstConfigs.plot(ax, color='red')

# Miscellaneous settings.
ax.set_xlim(1.00, 1.00e+11)
ax.set_ylim(1.00e+43, 1.00e+55)
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"$dE_{j}/dt$ [erg/s]")
ax.grid(True)
ax.legend()
plt.show()
