// File creates the python module in order to analyze GRB energy levels. Edit the file as you may see fit for other functions if needed.

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <cmath>

namespace py = pybind11;

const double pi = 3.14159265358979323846;
const double massOfProton = 1.67262192e-27;
const double c = 299792458;

// **************************************************************************

// Define parameters globally through a structure.
struct burstConfigurations {
    double initialEnergyEmitted; // Initial Energy Emitted by the GRB. Measured in erg.
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
    double sigma; // Precision Factor of Lorentz Factor phase transition model.
};

// Loads in values from the .yaml files.
burstConfigurations loadBurstConfigs(const py::dict &config_dict) {

    burstConfigurations burstConfigs; // Define the burstConfigurations structure via the "burstConfigs" alias.

    // Define variables in the .yaml files.
    burstConfigs.initialEnergyEmitted = config_dict["initialEnergyEmitted"].cast<double>();
    burstConfigs.initialLorentzFactor = config_dict["initialLorentzFactor"].cast<double>();
    burstConfigs.t_dec = config_dict["t_dec"].cast<double>();
    burstConfigs.t_NR = config_dict["t_NR"].cast<double>();
    burstConfigs.initialAngle = config_dict["initialAngle"].cast<double>();
    burstConfigs.initialMass = config_dict["initialMass"].cast<double>();
    burstConfigs.initialSolidAngle = config_dict["initialSolidAngle"].cast<double>();
    burstConfigs.r_NR = config_dict["r_NR"].cast<double>();
    burstConfigs.densityProportionalityConstant = config_dict["densityProportionalityConstant"].cast<double>();
    burstConfigs.powerLawIndex = config_dict["powerLawIndex"].cast<double>();
    burstConfigs.radiusProportionalityConstant = config_dict["radiusProportionalityConstant"].cast<double>();
    burstConfigs.surroundingMediumDensity = config_dict["surroundingMediumDensity"].cast<double>();
    burstConfigs.sigma = config_dict["sigma"].cast<double>();

    // Return the structure.
    return burstConfigs;

}

// **************************************************************************

// MAIN FUNCTIONS:

// **NOTE**
// Each functions will have an overloaded version of it that returns a double, in order to nest functions into other, more complex functions.
// The std::vector<double> functions are functions that will be exported into the actual pybind11 module, and will output an array to be graphed.
// Every important function will be exported given some alias at the bottom of this .cpp file to be used in python (i.e. lorentz_factor(t, burstConfigs)).

// **************************************************************************

double overloadedRadius(double t, const burstConfigurations &burstConfigs) {
    return burstConfigs.radiusProportionalityConstant * pow(t, 0.4);
}

// Function of radius from central engine wrt. time.
std::vector<double> radius(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {
        
        double y;
        y = overloadedRadius(t, burstConfigs);
        results.push_back(y);

    }

    return results;
}

// Overloaded version in order to nest the function into other functions.
double overloadedLorentzFactor(double t, const burstConfigurations &burstConfigs) {
    double xi = std::tanh(-burstConfigs.sigma * (t - burstConfigs.t_dec));
    return burstConfigs.initialLorentzFactor * (1 + xi) / 2.0 + burstConfigs.initialLorentzFactor * 
    (1 - xi) / (2.0 * std::sqrt(t / burstConfigs.t_dec));
}

std::vector<double> lorentzFactor(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {
        
        double y;
        y = overloadedLorentzFactor(t, burstConfigs);
        results.push_back(y);

    }
    return results;
}

// Overloaded version in order to nest the function into other functions.
double overloadedEjectionAngle(double t, const burstConfigurations &burstConfigs) {
    return burstConfigs.initialAngle + (burstConfigs.initialAngle / overloadedLorentzFactor(t, burstConfigs));
}

// Function of Ejection Angle of the GRB wrt. time.
std::vector<double> ejectionAngle(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {
        
        double y;
        y = overloadedEjectionAngle(t, burstConfigs);
        results.push_back(y);

    }
    return results;
}

// Power Law Density Model wrt. time.
double overloadedPowerLawDensityModel(double t, const burstConfigurations &burstConfigs) {
    return burstConfigs.densityProportionalityConstant * pow(overloadedRadius(t, burstConfigs), -1 * burstConfigs.powerLawIndex);
}

// Function of swept up mass density wrt. time.
double overloadedSweptUpMassDensity(double t, const burstConfigurations &burstConfigs) {
    return 4 * overloadedPowerLawDensityModel(t, burstConfigs) * overloadedLorentzFactor(t, burstConfigs);
}

// Calculates the derivative of the radius function in order to form the integral with respect to x.
double radiusDerivative(double x, const burstConfigurations &burstConfigs) {
        double stepSize = 1e-5;
        return (overloadedRadius(x + stepSize, burstConfigs) - overloadedRadius(x - stepSize, burstConfigs)) / (2 * stepSize);
}

double overloadedSweptUpMassPerSolidAngle(double t, const burstConfigurations &burstConfigs) {

    // Integral from limits 0.01 to t.
    int num_intervals = 100000;
    double integral = 0.0;
    double dx = (t - 0.01) / num_intervals;   

    for (int i = 0; i < num_intervals; ++i) {
        double x = 0.01 + i * dx;
        double x_next = x + dx;

        // Calculate integrand at x & x_next
        double integrand_At_x = overloadedLorentzFactor(x, burstConfigs) * overloadedSweptUpMassDensity(x, burstConfigs) * 
        pow(overloadedRadius(x, burstConfigs), 2) * radiusDerivative(x, burstConfigs);
        double integrand_At_x_next = overloadedLorentzFactor(x_next, burstConfigs) * 
        overloadedSweptUpMassDensity(x_next, burstConfigs) * pow(overloadedRadius(x_next, burstConfigs), 2) * 
        radiusDerivative(x_next, burstConfigs);

        // Implement trapezoidal rule.
        integral += 0.5 * (integrand_At_x + integrand_At_x_next) * dx;
    }

    return integral;
}

// Calculate swept up mass per solid angle of the GRB.
std::vector<double> sweptUpMassPerSolidAngle(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y = overloadedSweptUpMassDensity(t, burstConfigs);
        results.push_back(y);

    }

    return results;
}

double overloadedEjectaMassPerSolidAngle(double t, const burstConfigurations &burstConfigs) {
    return overloadedSweptUpMassPerSolidAngle(t, burstConfigs) + (burstConfigs.initialMass / burstConfigs.initialSolidAngle);
}

// Calculate ejecta mass per solid angle by referencing the overloaded function above.
std::vector<double> ejectaMassPerSolidAngle(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y;
        y = overloadedEjectaMassPerSolidAngle(t, burstConfigs);
        results.push_back(y);

    }
    return results;
}

// **************************************************************************

// OVERLOADED VERSIONS OF JET ENERGY FUNCTION

// **************************************************************************

// Function of alpha, detailed in the paper.
double alpha(double t, const burstConfigurations &burstConfigs) {
    return pow((pow(
    (overloadedEjectaMassPerSolidAngle(t, burstConfigs) / overloadedSweptUpMassPerSolidAngle(t, burstConfigs)), 0.25)) * 
    (pow((overloadedLorentzFactor(t, burstConfigs) / burstConfigs.initialLorentzFactor), 0.5)) *
    (pow((burstConfigs.initialAngle / overloadedEjectionAngle(t, burstConfigs)), 0.5)), 0.5);
}

// Function of beta, detailed in the paper.
double beta(double t, const burstConfigurations &burstConfigs) {
    return 1 / ((pow(burstConfigs.initialLorentzFactor, 0.5)) * (pow((burstConfigs.r_NR / overloadedRadius(t, burstConfigs)), 1) * 
    pow((burstConfigs.initialEnergyEmitted / (burstConfigs.initialEnergyEmitted + 1.5e+50)), 1) * 
    (pow((overloadedEjectaMassPerSolidAngle(t, burstConfigs) / overloadedSweptUpMassPerSolidAngle(t, burstConfigs)), 0.25))));
} 

// Evaluate the derivative of the overloaded jet energy function with respect to time, BEFORE the non-relativistic stage.
double overloadedJetEnergy1(double t, const burstConfigurations &burstConfigs) {
   return 1.00e+45 + (burstConfigs.initialEnergyEmitted - 1.00e+45) * (1 / (t * (1 + pow(t / burstConfigs.t_dec, alpha(t, burstConfigs)))));
}

double overloadedJetEnergy2(double t, const burstConfigurations &burstConfigs) {
    return overloadedJetEnergy1(burstConfigs.t_NR, burstConfigs) * std::exp(-1 * beta(t, burstConfigs) * (t - burstConfigs.t_NR));
} 

// **************************************************************************

// FIRST MAIN FUNCTIONS OF THE PAPER: JET ENERGY LEVEL WITH RESPECT TO TIME (s) BEFORE [= jetEnergy1] AND AFTER [= jetEnergy2] t = t_NR.

// **************************************************************************

std::vector<double> jetEnergy1(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y;
        y = (overloadedJetEnergy1(t, burstConfigs));
        results.push_back(y);

    }
    return results;
}

std::vector<double> jetEnergy2(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y;
        y = (overloadedJetEnergy2(t, burstConfigs));
        results.push_back(y);

    }
    return results;
}

// **************************************************************************

// ENERGY PER SOLID ANGLE TESTS

// **************************************************************************

// Calculate overloaded version of Energy Per Solid Angle functions
double overloadedEnergyPerSolidAngle1(double t, const burstConfigurations &burstConfigs) {
   return burstConfigs.initialEnergyEmitted;
}

double overloadedEnergyPerSolidAngle2(double t, const burstConfigurations &burstConfigs) {
    return (overloadedEnergyPerSolidAngle1(burstConfigs.t_NR, burstConfigs) * std::exp(-beta(t, burstConfigs) * pow((t - burstConfigs.t_NR), 1.3))) *
    (1 + pow((t - burstConfigs.t_NR) / burstConfigs.t_dec, pow(alpha(t, burstConfigs), 5)));
} 

// **************************************************************************

// SECOND MAIN FUNCTIONS OF THE PAPER: ENERGY PER SOLID ANGLE WITH RESPECT TO TIME (s).

// ************************************************************************** 

std::vector<double> energyPerSolidAngle1(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y;
        y = overloadedEnergyPerSolidAngle1(t, burstConfigs);
        results.push_back(y);

    }
    return results;
}

std::vector<double> energyPerSolidAngle2(const std::vector<double> &t_values, const burstConfigurations &burstConfigs) {

    std::vector<double> results;
    results.reserve(t_values.size());

    for (double t : t_values) {

        double y;
        y = overloadedEnergyPerSolidAngle2(t, burstConfigs);
        results.push_back(y);

    }
    return results;
}

// **************************************************************************

// EXTERNAL MODELS

// **************************************************************************

std::vector<double> gaussianModel1(const std::vector<double> &theta_values) {

    std::vector<double> results;
    results.reserve(theta_values.size());
    double thetaC = 3;

    for (double theta : theta_values) {

        double y;
        y = (1e+53 / (4 * pi)) * std::exp((-1 * pow((theta / thetaC), 2) / 2));
        results.push_back(y);

    }
    return results;
}

// INDIVIDUAL RATIO TESTS FOR GENERAL JET ENERGY ANALYSIS AND TROUBLESHOOTING.
// ****************************************************************************


// DEBUG:
// - JET ENERGY FUNCTION ** DONE ** 
// - ENERGY PER SOLID ANGLE FUNCTION
// - RADIUS RATIO ** DONE ** 

// Package the functions into the JET_ENERGY_FUNCTIONS module.
PYBIND11_MODULE(JET_ENERGY_FUNCTIONS, m) {

    // Import the burst configuration variables.
    py::class_<burstConfigurations>(m, "burstConfigurations")

        .def(py::init<>())
        .def_readwrite("initialEnergyEmitted", &burstConfigurations::initialEnergyEmitted)
        .def_readwrite("initialLorentzFactor", &burstConfigurations::initialLorentzFactor)
        .def_readwrite("t_dec", &burstConfigurations::t_dec)
        .def_readwrite("t_NR", &burstConfigurations::t_NR)
        .def_readwrite("initialAngle", &burstConfigurations::initialAngle)
        .def_readwrite("initialMass", &burstConfigurations::initialMass)
        .def_readwrite("initialSolidAngle", &burstConfigurations::initialSolidAngle)
        .def_readwrite("r_NR", &burstConfigurations::r_NR)
        .def_readwrite("densityProportionalityConstant", &burstConfigurations::densityProportionalityConstant)
        .def_readwrite("powerLawIndex", &burstConfigurations::powerLawIndex)
        .def_readwrite("radiusProportionalityConstant", &burstConfigurations::radiusProportionalityConstant)
        .def_readwrite("surroundingMediumDensity", &burstConfigurations::surroundingMediumDensity)
        .def_readwrite("alpha", &burstConfigurations::sigma);

    // Create the functions using pybind11.
    m.def("load_burst_configs", &loadBurstConfigs, "Loads the GRB configurations from the .yaml files.");
    m.def("radius", &radius, "A function of radius from the central engine of the GRB in centimeters (cm).");
    m.def("lorentz_factor", &lorentzFactor, "A function that finds the approximate Lorentz factor of a GRB.");
    m.def("ejection_angle", &ejectionAngle, "A function that finds the approximate ejection angle of the GRB in radians.");
    m.def("mass_per_solid_angle", &sweptUpMassPerSolidAngle, "A function that calculates the mass swept up per solid angle over the GRB's lifespan.");
    m.def("ejecta_mass_per_solid_angle", &ejectaMassPerSolidAngle, "A function that calculates the total ejecta mass of the GRB with respect to time, in seconds (s).");

    m.def("JET_ENERGY1", &jetEnergy1, "A function that represents the approximate energy of a GRB jet before t = t_NR.");
    m.def("JET_ENERGY2", &jetEnergy2, "A function that represents the approximate energy of a GRB jet after t = t_NR.");

    m.def("ENERGY_PER_SOLID_ANGLE1", &energyPerSolidAngle1, "A function that represents the approximate energy of a GRB jet before t = t_NR.");
    m.def("ENERGY_PER_SOLID_ANGLE2", &energyPerSolidAngle2, "A function that represents the approximate energy of a GRB jet after t = t_NR.");

    m.def("gaussian_model", &gaussianModel1);
}