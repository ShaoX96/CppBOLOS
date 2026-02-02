/* Master code of CppBOLOS.
*  An example of H2/O2/He mixture under constant reduced E-field.
*
* This file is part of CppBOLOS. See license and copyright information at
* https://github.com/ShaoX96/CppBOLOS
*/

#include <iostream>
#include "grid.h"
#include "solver.h"
#include <vector>
#include <functional>
#include <sstream>
#include "Logger.h"
#include "utilities.h"

using namespace CppBOLOS;

int main() {

    currentLogLevel = LOG_DEBUG; // LOG_NONE, LOG_WARNING

    // Read cross section data
    std::string filename = "../data/bolsigdb_H2O2He.dat";
    std::vector<std::string> ss = clean_file(filename);
    std::vector<Collision> collisions = parse(ss);

    // Create a BoltzmannSolver instance
    BoltzmannSolver bsolver;

    // Set up grid
    std::unique_ptr<LinearGrid> linGrid(new LinearGrid(0.0, 80, 200));
    bsolver.set_grid(std::move(linGrid));
    // Alternatively, user can set up grid simply by:
    // bsolver.set_grid("QuadraticGrid", 0, 60, 180);

    // Load collision data
    bsolver.load_collisions(collisions);
    LOG_INFO("\nA total of " + std::to_string(bsolver.number_of_targets()) +
             " targets have been loaded:\n" + bsolver.targetNames());

    // Set density, T_gas, E/N, and Initialization
    double Tgas = 300;    // Gas teperature [K]
    double EN = 300;      // Reduced electirc field [Td]
    bsolver.set_density("H2:0.1667, O2:0.0833, HE:0.75");
    bsolver.set_kT (Tgas);
    bsolver.set_EN(EN);
    bsolver.init();

    // Start from Maxwell EEDF
    Eigen::VectorXd f0 = bsolver.maxwell(EN/50); // a rough of estimation of mean electron energy

    // Converge
    Eigen::VectorXd f1;
    try {
        f1 = bsolver.converge(f0, 200, 1e-4);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    }

    // Check f1 values
    vecPrint(f1);

    // Print some useful parameter
   std::cout << "Rate O2 -> O2^+: " << bsolver.rate(f1,"O2->O2^+") << std::endl;

    std::cout << "reduced mobility: " << bsolver.mobility() << std::endl;

    std::cout << "reduced diffusion: " << bsolver.diffusivity() << std::endl;

    double mean_energy = bsolver.mean_energy(f1);
    std::cout << "mean energy: " << mean_energy << "eV (" << bsolver.get_Te() << "K)" << std::endl;

    std::cout << "EN: " << bsolver.get_EN() << std::endl;

    std::cout << "Power/N [eV m^3/s]: " << bsolver.elec_power() << std::endl;

    // Second pass
    Grid& oldGrid = bsolver.get_grid(); // set old grid
    bsolver.set_grid("QuadraticGrid", 0.0, mean_energy*15, 200); // set new grid
    Grid& newGrid = bsolver.get_grid(); // get new grid
    bsolver.init();

    auto f2 = newGrid.interpolate(f1, oldGrid);  // interpolate old solution & grid on new grid
    f2 = bsolver.converge(f2, 200, 1e-5);
    vecPrint(f2);
    mean_energy = bsolver.mean_energy(f2);
    std::cout << "mean energy: " << mean_energy << "eV (" << bsolver.get_Te() << "K)" << std::endl;

    return 0;
}


