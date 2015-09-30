
// Set the number of atoms in the box
const int n_atoms = 1000;

// Set the number of Monte Carlo moves to perform
//const int eq_moves  = 10000000;
const int eq_moves  = 2500;
const int num_moves =  2000000;

// Set the size of the box (in Angstroms)
const double box_size[3] = { 100.0, 100.0, 100.0 };

// The maximum amount that the atom can be translated by
//const double max_translate = 0.4e-2;  // angstroms
//const double max_translate = 0.01;  // angstroms
const double max_translate = 0.0023;  // angstroms

// Simulation temperature
const double temperature = 273.15;   // kelvin

// calculate kT
const double k_boltz = 1.987206504191549E-003;  // kcal mol-1 K-1
const double kT = k_boltz * temperature;

// Give the Lennard Jones parameters for the atoms
// (these are the OPLS parameters for Krypton)
double sigmaP = 1.0000;     // angstroms
double sigmaI = 0.0001;   // angstroms ICE
const double epsilon = 0.317;   // kcal mol-1
const int nsites = 5; //Nucleation Sites
const double dt = 0.0005;
const double D  = 0.2235;
const double Ddt = 2.0 * D * dt;
//const double IceGrowth = 0.1e-3;
const double IceGrowth = 1.0 *dt;
