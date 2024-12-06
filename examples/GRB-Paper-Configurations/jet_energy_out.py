import numpy as np
import matplotlib.pyplot as plt
import yaml
import JET_ENERGY_FUNCTIONS as JEF  # Imports the compiled C++ module containing all the functions needed for analysis.

class GRB_jet_energy:
    def __init__(self, name, config_path):
        self.name = name
        self.config_path = config_path
        self.parameters = self._load_yaml(f"{config_path}/Parameters-{name}.yaml")
        self.calculated_parameters = self._load_yaml(f"{config_path}/CalculatedParameters-{name}.yaml")
        self.burst_configs = JEF.load_burst_configs({**self.parameters, **self.calculated_parameters})
        self.t_values1 = np.linspace(1.00, self.burst_configs.t_NR, 500).tolist()
        self.t_values2 = np.linspace(self.burst_configs.t_NR, 1.00e+11, 500).tolist()

    def _load_yaml(self, file_path):
        with open(file_path) as file:
            return yaml.safe_load(file)

    def plot(self, ax, color):
        # Plot jet energy for time ranges.
        ax.loglog(self.t_values1, JEF.JET_ENERGY1(self.t_values1, self.burst_configs), color=color, label=f"{self.name}")
        ax.loglog(self.t_values2, JEF.JET_ENERGY2(self.t_values2, self.burst_configs), color=color)
        # Plot t = t_NR line.
        ax.axvline(x=self.burst_configs.t_NR, color=color, linestyle='--', label=fr"$t = t_{{NR, {self.name[-1]}}}$")

# Define Plotted GRBs.
grb1 = GRB_jet_energy(name="GRB1", config_path="C:/repos/GRB-Jet-Energy-Model/examples/GRB-Paper-Configurations/config/GRB1-Configurations")
grb2 = GRB_jet_energy(name="GRB2", config_path="C:/repos/GRB-Jet-Energy-Model/examples/GRB-Paper-Configurations/config/GRB2-Configurations")
grb3 = GRB_jet_energy(name="GRB3", config_path="C:/repos/GRB-Jet-Energy-Model/examples/GRB-Paper-Configurations/config/GRB3-Configurations")

# Plot Each GRB.
fig, ax = plt.subplots(figsize=(6, 7))
grb1.plot(ax, color='red')
grb2.plot(ax, color='blue')
grb3.plot(ax, color='orange')

# Miscellaneous settings.
ax.set_xlim(1.00, 1.00e+11)
ax.set_ylim(1.00e+43, 1.00e+55)
ax.set_xlabel("Time [s]")
ax.set_ylabel(r"$dE_{j}/dt$ [erg/s]")
ax.grid(True)
ax.legend()
plt.show()
