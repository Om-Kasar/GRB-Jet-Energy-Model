# RUN THE FOLLOWING FILE AFTER EDITING THE "Parameters.yaml" FILE IN THE DIRECTORY.

import numpy as np
import math
import yaml

c = 2.99792458e+10 # Speed of light, in cm/s.
protonMass = 1.67262192e-24 # Mass of a proton, in g.

with open("C:/repos/GRB-Jet-Energy-Model/config/Parameters.yaml", "r") as file: # Put your path to the Parameters.yaml file.
    config = yaml.safe_load(file)

t_dec = float(config["t_dec"])
density_Proportionality_Constant = float(config["densityProportionalityConstant"])
power_Law_Index = float(config["powerLawIndex"])
radius_Proportionality_Constant = float(config["radiusProportionalityConstant"])
surrounding_Medium_Density = float(config["surroundingMediumDensity"])
initial_Energy_Emitted = float(config["initialEnergyEmitted"])
pi = math.pi

# Days to Seconds function, for t_NR
def DaysToSeconds(dayvalue):
    return dayvalue * 24 * 60 * 60

# Solar Mass to G, if needed.
def SolarMasstoG(solarMass):
    return solarMass * 1.989e33

def RadiansToDegrees(radianMeasure):
    return radianMeasure * (180 / pi)

# Calculate Parameters
# ************************

initialAngle = float(input("Give the initial ejection angle of the GRB, in radians: ")) # Initial Angle of GRB, in radians. Input degree/radian value for Initial Angle, use the DegreesToRadians function if needed.
initialIsotropicEquivalentEnergy = (2 * initial_Energy_Emitted) / (initialAngle ** 2)
initialLorentzFactor = (initialIsotropicEquivalentEnergy / (4 * pi * protonMass * surrounding_Medium_Density *  (c ** 5) * (t_dec**3)))**(1/8)
initialMass = initial_Energy_Emitted / (initialLorentzFactor * c**2)
initialSolidAngle = 2 * pi * (1 - math.cos(initialAngle))
t_NR = DaysToSeconds((1100 * (initialIsotropicEquivalentEnergy / (1.00e+53 * surrounding_Medium_Density))) ** (1/3))
r_NR = radius_Proportionality_Constant * (t_NR ** 2/5)
sigma = 5 / t_dec # Smoothing Factor

# Update values to the dynamically created CalculatedParameters.yaml file.

New_Parameters = {
    "initialAngle": initialAngle,
    "initialIsotropicEquivalentEnergy": initialIsotropicEquivalentEnergy,
    "initialLorentzFactor": initialLorentzFactor,
    "initialMass": initialMass,
    "initialSolidAngle": initialSolidAngle,
    "t_NR": t_NR,
    "r_NR": r_NR,
    "sigma": sigma,
    "pi": pi,
    "protonMass": protonMass,
    "c": c
}

with open("C:/repos/GRB-Jet-Energy-Model/config/CalculatedParameters.yaml", "w") as file1:

    file1.write("# Dynamically created parameters for calculations. \n")

    yaml.dump(New_Parameters, file1)

print("Calculations complete. Results saved to CalculatedParameters.yaml.")
print("**DELETE THE 'CalculatedParameters.yaml' FILE WHEN TESTING & EDITING OTHER GRBS, OR IF ANY PROBLEMS OCCUR.**")