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
t_values = np.linspace(1.00, 1.00e+10, 500).tolist()

# Plot the values.
plt.figure(figsize=(10, 6))
plt.loglog(t_values1, JEF.JET_ENERGY1(t_values1, burstConfigs), color = 'blue')
plt.loglog(t_values2, JEF.JET_ENERGY2(t_values2, burstConfigs), color = 'blue')
plt.axvline(x = burstConfigs.t_NR, color = 'red', linestyle = '--', label = "$t = t_{NR}$")
plt.ylim(1.00e+43, 1.00e+51)

plt.xlabel("Time [s]")
plt.ylabel("Energy [erg]")
plt.grid(True)
plt.legend()
plt.show()