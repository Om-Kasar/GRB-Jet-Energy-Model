## Other Usages

In order to display multiple graphs (similar to that of the referenced research paper), you need to create however many configuration folders needed for comparisons. Refer to the research paper example folder if needed.

1. Make seperate config folder containing files for each GRB being compared. Utilize a "parameters.py" file in order to make the "CalculatedParameters.py" file for each GRB being analyzed. Take a look at example 2 if needed.

2. After creating each GRB's configurations in their respective folders, create multiple instances of the "burstConfigs" for each GRB, using the C++ system's "jef.load_burst_configs()" function.

```python
import jet_energy_functions as jef

# Load GRB configurations for each GRB (modify paths as needed)
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
```

3. Create a figure to analyze the distribution of jet energy using matplotlib, and plot it accordingly.

```python
# Create a figure for Jet Energy plot
fig1, ax1 = plt.subplots(figsize=(6, 7))

# Create instances of classes for Jet Energy
jet_energy_calculator1 = jef.jet_energy(burstConfigs1)
jet_energy_calculator2 = jef.jet_energy(burstConfigs2)
jet_energy_calculator3 = jef.jet_energy(burstConfigs3)

jet_energy_calculator1.plot(ax1, color='red', label=1)
jet_energy_calculator2.plot(ax1, color='blue', label=2)
jet_energy_calculator3.plot(ax1, color='orange', label=3)

# Feel free to edit the graph's configurations as needed.
```

4. Create a seperate figure that analyzes the energy per solid angle of the GRB. Make edits to the graph's x and y-axis accordingly.

```python
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
```

5. Finally, run the program. Images like this should be outputted for both figures:

![image](C:/GRB paper images/Figure_3.png)
![image](C:/GRB paper images/Figure_4.png)