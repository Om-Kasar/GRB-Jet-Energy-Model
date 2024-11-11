/*****************************************************************
Authors: Om Kasar & Annika Thomas
Description: This program defines the functions needed to analyze energy distributions of GRB jets. 
Link to paper: **LINK TO PAPER HERE**
*****************************************************************/
#include <iostream>
#include <cmath>
#include <functional>
#include <yaml-cpp/yaml.h>
#include <pybind11/pybind11.h>
using std::tanh;

namespace py = pybind11;

const double pi = 3.14159265358979323846;
const double massOfProton = 1.67262192e-27;
const double c = 299792458;

// IMPORTANT: The 't' parameter in most of the functions represents time, in seconds (s).

// Define parameters globally.
struct Params {
    double initialEnergyEmitted; // Initial Energy Emitted by the GRB. Measured in erg.
    double totalEnergyEmitted; // Energy Emitted from a given closed interval [0, t].
    double initialLorentzFactor; // Initial Lorentz Factor of the GRB. Unitless value.
    double t_dec; // Deceleration Time of the GRB. Measured in seconds (s).
    double t_NR; // Time in which the GRB becomes non-relativistic. Measured in seconds (s).
    double initialAngle; // Initial Angle (theta) of the GRB. Measured in radians (rad).
    double initialMass; // Initial Mass of the GRB. Needed in order to evaluate the function of total swept up mass. Measured in grams (g).
    double initialSolidAngle; // Initial Solid Angle of the GRB. Used for evaluating the intial mass over the initial solid angle. Measured in steradians (str).
    double r_NR; // The radius length in which the GRB transitions into the non-relativistic stage.
    double densityProportionalityConstant; // Proportionality constant to evaluate the mass swept up density function.
    double powerLawIndex; // Power Law Index parameter in order to model the power law density model.
    double radiusProportionalityConstant; // Proportionality constant to evaluate the function of radius.
    double surroundingMediumDensity; // Surrounding Medium Density of the interstellar environment.
    double initialIsotropicEquivalentEnergy; // Initial Isotropic Equivalent Energy of the GRB.
    double alpha; // Precision Factor of Lorentz Factor phase transition model.
};

Params globalParams;

// This work references these global parameters in the structure (containing values from the .yaml files) multiple times in most of the functions, as suggested by the 'globalParams.' 
// prefix in front of each parameter.

void loadGlobalParameters() {
    // Load the two .yaml files located in the config folder.
    YAML::Node Parameters = YAML::LoadFile("C:/POLYGENCEPROJECT/config/Parameters.yaml"); // Put your path to the Parameters.yaml file here.
    YAML::Node CalculatedParameters = YAML::LoadFile("C:/POLYGENCEPROJECT/config/CalculatedParameters.yaml"); // Put your path to the CalculatedParameters.yaml file here.

    // Call values from the Parameters.yaml file, located in the config folder.
    globalParams.t_dec = Parameters["Params"]["t_dec"].as<double>();
    globalParams.densityProportionalityConstant = Parameters["Params"]["densityProportionalityConstant"].as<double>();
    globalParams.powerLawIndex = Parameters["Params"]["powerLawIndex"].as<double>();
    globalParams.radiusProportionalityConstant = Parameters["Params"]["radiusProportionalityConstant"].as<double>();
    globalParams.surroundingMediumDensity = Parameters["Parameters"]["surroundingMediumDensity"].as<double>();
    globalParams.initialEnergyEmitted = Parameters["Params"]["initialEnergyEmitted"].as<double>();

    // Call values from the CalculatedParameters.yaml file, located in the config folder.
    globalParams.initialLorentzFactor = CalculatedParameters["CalculatedParams"]["initialLorentzFactor"].as<double>();
    globalParams.initialAngle = CalculatedParameters["CalculatedParams"]["initialAngle"].as<double>();
    globalParams.initialIsotropicEquivalentEnergy = CalculatedParameters["CalculatedParams"]["initialIsotropicEquivalentEnergy"].as<double>();
    globalParams.initialMass = CalculatedParameters["CalculatedParams"]["initialMass"].as<double>();
    globalParams.initialSolidAngle = CalculatedParameters["CalculatedParams"]["initialSolidAngle"].as<double>();
    globalParams.r_NR = CalculatedParameters["CalculatedParams"]["r_NR"].as<double>();
    globalParams.t_NR = CalculatedParameters["CalculatedParams"]["r_NR"].as<double>();

    // Variable that declares the level of accuracy of the hyperbolic tangent Lorentz factor function.
    globalParams.alpha = 0.4; // As the value approaches zero, the approximation of the Lorentz factor function becomes more and more accurate to the Phase Transition model.

}

// Function of radius from central engine wrt. time.
double radius(double t) {
    return globalParams.radiusProportionalityConstant * pow(t, 2/5);
}

// Function of Lorentz Factor of the GRB wrt. time, using an approximation of the Phase Transition model via hyperbolic tangent.
double lorentzFactor(double t) {
    return globalParams.initialLorentzFactor * pow(globalParams.t_dec / t, 0.5) *
           (1 + tanh(globalParams.alpha * (t - globalParams.t_dec))) / 2 +
           globalParams.initialLorentzFactor * (1 - tanh(globalParams.alpha * (t - globalParams.t_dec))) / 2;
}

// Function of Ejection Angle of the GRB wrt. time.
double ejectionAngle(double t) {
    return globalParams.initialAngle + (globalParams.initialAngle / lorentzFactor(t));
}

// Power Law Density Model wrt. time.
double powerLawDensityModel(double t) {
    return globalParams.densityProportionalityConstant * (pow(radius(t), -1 * globalParams.powerLawIndex));
}

// Density of Swept Up Mass wrt. time.
double sweptUpMassDensity(double t) {
    return 4 * powerLawDensityModel(t) * lorentzFactor(t);
}

// Implimenting differentiation.
double differentiate(double (*func)(double), double x) {
    double stepSize = 1e-6;
    return (func(x + stepSize) - func(x - stepSize)) / (2 * stepSize);
}

// Implementing integration.
double integrate(std::function<double(double)> func, 
double lowerLimit, double upperLimit, int numOfSubintervals) {

    double stepSize = (upperLimit - lowerLimit) / numOfSubintervals;
    double sum = 0.5 * (func(lowerLimit) + func(upperLimit));

    for (int i = 1; i < numOfSubintervals; ++i) {
        double x = lowerLimit + i * stepSize;
        sum += func(x);
    }

    return sum * stepSize;

}

// Define the integrand.
double INTEGRAND(double x) {
    return lorentzFactor(x) * sweptUpMassDensity(x) * pow(radius(x), 2) * differentiate(radius, x);
}

// Mass Per Solid Angle wrt. time.
double sweptUpMassPerSolidAngle(double t) {
        // Finds the integrand using a dummy variable 'x' in order to integrate with respect to dx.
        // The dummy variable 'x' is used in order to plug in some value t into the integrated function afterwards.
        
        double result = integrate(INTEGRAND, 0, t, 10000);
    
    return result;
}

// Total Swept Up Mass per Solid Angle of GRB wrt. time.
double totalMassPerSolidAngle(double t) {
    return sweptUpMassPerSolidAngle(t) + (globalParams.initialMass / globalParams.initialSolidAngle);
     // Adding ratio of initial mass & solid angle for the units to comply with units.
}

// Function of Solid Angle wrt. time.
double solidAngle(double t) {
    return 2 * pi * (1-cos(ejectionAngle(t)));
}

// MAIN FUNCTION OF THE PAPER: Jet Energy wrt. time.
double jetEnergy(double t) {
    try {
        if (0 < t < globalParams.t_NR) {
            return globalParams.initialEnergyEmitted * (totalMassPerSolidAngle(t) / sweptUpMassPerSolidAngle(t)) * (pow((lorentzFactor(t) / globalParams.initialLorentzFactor), 1.5)) *
            (pow((globalParams.initialEnergyEmitted / ejectionAngle(t)), 0.5));
        } else if (t >= globalParams.t_NR) {
            return globalParams.initialEnergyEmitted * (totalMassPerSolidAngle(t) / sweptUpMassPerSolidAngle(t)) * (pow((globalParams.r_NR / radius(t)), 1.5)) * 
            (pow((globalParams.initialEnergyEmitted / 1e+54), 0.5));
        }
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }
}

//  Calculations for Energy Per Solid Angle Analysis:
// ****************************************************

// Derivative of Jet Energy with respect to time.
double jetEnergyDerivative(double t) {
    double result = differentiate(jetEnergy, t);
    return result;
}

// Derivative of Solid Angle formula, with respect to time. 
double solidAngleDerivative(double t) {
    double result = differentiate(solidAngle, t);
    return result; // Gives final result.

}

// Gives energy per solid angle function.
double energyPerSolidAngle(double t) {
    return jetEnergyDerivative(t) * (1 / (solidAngleDerivative(t) * differentiate(ejectionAngle, t))); // Energy Per Solid Angle [dE/dt * dt/d(Omega)].
}

// Compile pybind11 functions.

PYBIND11_MODULE(computationalfunctions, m) {

    m.doc() = "This library imports all the functions needed in order to approximate the energy distribution of GRB jets.";

    m.def("lorentz_factor", &lorentzFactor, "A function that finds the approximate Lorentz factor of a GRB.");
    m.def("ejection_angle", &ejectionAngle, "A function that finds the approximate ejection angle of the GRB.");
    m.def("radius", &radius, "A function of radius from the central engine of the GRB.");
    m.def("power_law_model", &powerLawDensityModel, "A standard power law model of a GRB.");
    m.def("swept_up_mass_density", &sweptUpMassDensity, "A function illustrating the density of swept up mass.");
    m.def("mass_per_solid_angle", &sweptUpMassPerSolidAngle, "A function that calculates the mass swept up per solid angle over the GRB's lifespan.");
    m.def("ejecta_mass_per_solid_angle", &totalMassPerSolidAngle, "A function that resperents the total ejecta mass.");
    m.def("jet_energy", &jetEnergy, "A function that represents the approximate energy of a GRB jet.");
    m.def("solid_angle", &solidAngle, "A function of Solid Angle.");
    m.def("jet_energy_derivative", &jetEnergyDerivative, "The derivative of the jet energy function, used for energy per solid angle analysis.");
    m.def("solid_angle_derivative", &solidAngleDerivative, "The derivative of the solid angle function, used for energy per solid angle analysis.");
    m.def("energy_per_solid_angle", &energyPerSolidAngle, "A function of energy per solid angle analysis.");

}