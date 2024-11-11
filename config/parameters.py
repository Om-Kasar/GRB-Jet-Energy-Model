import numpy as np
import math
import yaml

with open("C:/POLYGENCEPROJECT/config/Parameters.yaml", "r") as file: # Put your path to the Parameters.yaml file.
    config = yaml.safe_load(file)

c = float(config["Params"]["c"])
proton_Mass = float(config["Params"]["protonMass"])
pi = float(config["Params"]["pi"])
t_dec = float(config["Params"]["t_dec"])
density_Proportionality_Constant = float(config["Params"]["densityProportionalityConstant"])
power_Law_Index = float(config["Params"]["powerLawIndex"])
radius_Proportionality_Constant = float(config["Params"]["radiusProportionalityConstant"])
surrounding_Medium_Density = float(config["Params"]["surroundingMediumDensity"])
initial_Energy_Emitted = float(config["Params"]["initialEnergyEmitted"])

# Days to Seconds function, for t_NR
def DaysToSeconds(dayvalue):
    return dayvalue * 24 * 60 * 60

# Solar Mass to G, if needed.
def SolarMasstoG(solarMass):
    return solarMass * 1.989e33

# Degrees to Radians Function
def DegreesToRadians(degreeMeasure):
    return degreeMeasure * (np.pi / 180)

# Calculate Parameters
initialAngle = DegreesToRadians(41.5073) # Initial Angle of GRB. Input degree/radian value for Initial Angle, use the DegreesToRadians function if needed.
initialIsotropicEquivalentEnergy = (2 * initial_Energy_Emitted) / (initialAngle) # Calculates the initial Isotropic Equivalent Energy of the GRB.
initialLorentzFactor = (initialIsotropicEquivalentEnergy / (4 * pi * proton_Mass * surrounding_Medium_Density *  (c ** 5) * (t_dec**3)))**(1/8) # Calculates initial Lorentz Factor.
initialMass = initial_Energy_Emitted / (initialLorentzFactor * c**2) # Calculates initial Lorentz Factor. Refer to the paper for details.
initialSolidAngle = 2 * pi * (1 - math.cos(initialAngle)) # Calculates initial Solid Angle
t_NR = DaysToSeconds(1100 * (initialIsotropicEquivalentEnergy / 1e+53 * surrounding_Medium_Density)) # Calculates t_NR, the time it takes for the GRB to become non-relativistic.
r_NR = radius_Proportionality_Constant * (t_NR ** 2/5) # Calculates the radius from the central engine at which the GRB becomes non-relativistic [ = r(t_NR)].

# Update values to the new CalculatedParameters.yaml file.

New_Parameters = {
    "CalculatedParams": {
        "initialAngle": initialAngle,
        "initialIsotropicEquivalentEnergy": initialIsotropicEquivalentEnergy,
        "initialLorentzFactor": initialLorentzFactor,
        "initialMass": initialMass,
        "initialSolidAngle": initialSolidAngle,
        "t_NR": t_NR,
        "r_NR": r_NR
    }
}


with open("C:/POLYGENCEPROJECT/config/CalculatedParameters.yaml", "w") as file1:

    file1.write("# Dynamically created parameters for calculations. \n")

    yaml.dump(New_Parameters, file1)


print("Calculations complete. Results saved to CalculatedParameters.yaml.")
print("**DELETE THE 'CalculatedParameters.yaml' FILE WHEN TESTING & EDITING OTHER GRBS, OR IF ANY PROBLEMS OCCUR.**")