#ifndef SOLVER_H
#define SOLVER_H

#include "grid.h"
#include "process.h"
#include "target.h"
#include "parser.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <variant>
#include <memory> // Needed for std::unique_ptr
#include <stdexcept> // For std::runtime_error
#include <cmath>
#include "Eigen/Sparse"
#include "Eigen/Dense"

namespace CppBOLOS {

// Constants
constexpr double ELEMENTARY_CHARGE = 1.602176634e-19;  // Coulombs
constexpr double ELECTRON_MASS = 9.10938356e-31;  // kg
const double GAMMA = std::sqrt(2 * ELEMENTARY_CHARGE / ELECTRON_MASS);
constexpr double TOWNSEND = 1e-21;
constexpr double KB = 1.380649e-23;  // J/K
constexpr double ELECTRONVOLT = 1.602176634e-19;  // Joules
constexpr double EV2K = 1.16045052e4;

// Equivalent to Python's "ConvergenceError" exception class
class ConvergenceError : public std::runtime_error {
public:
    ConvergenceError(const std::string& message)
            : std::runtime_error("CppBOLOS::ConvergenceError: " + message) {}
};

// TODO: You may need to define or import additional types or libraries.

class BoltzmannSolver {
public:
//    BoltzmannSolver(Grid* grid);  // Constructor takes a Grid object
    BoltzmannSolver() {} // default constructor

    std::string targetNames() const {
        std::string targetsInfo;
        for (const auto& pair : m_targetsMap) {
            targetsInfo += pair.first + " ";
        }
        return targetsInfo;
    }

    size_t number_of_targets() const { return m_targetsMap.size(); }

    // TODO: Add appropriate constructors or methods if you need to set up
    // the BoltzmannSolver object in a specific way.

    void set_grid(std::unique_ptr<Grid> grid);
    void set_grid(const std::string& type, double start, double end, int N_cells);
    void set_grid() {
        std::cout << "[Alert] Using default grid settings: \"LinearGrid\", (0.0, 60.0, 200) \n";
        set_grid("LinearGrid", 0.0, 60.0, 200);
    }

     Grid& get_grid() const {
        if (!m_grid) {
            throw std::runtime_error("Grid has not been set");
        }
        return *m_grid;
    }

    Eigen::VectorXd get_gridCell() const { return this->cenergy; }

    // TODO: You may need to add a destructor if there's any cleanup required
    // when a BoltzmannSolver object is destroyed.
    double get_kT() const { return kT * ELECTRONVOLT / KB; }
    void set_kT(double value) { kT = value * KB / ELECTRONVOLT; }
    double get_EN() const { return EN / TOWNSEND; }
    void set_EN(double value) { EN = value * TOWNSEND; }

    void set_density(const std::string& species, const double& density);
    void set_density(const std::string& ss);
    void set_density(const std::map<std::string, double>& xMap);
    std::vector<Process> load_collisions(const std::vector<Collision>& collisions);
    // Process add_process(std::map<std::string, std::variant<std::string, double, std::vector<std::pair<double, double>>>> kwargs);
    Process add_process(const std::string& target, const std::string& kind,
                        const std::vector<std::vector<double>>& data,
                        const double &mass_ratio, const std::string &product,
                        const double &threshold, const double &weight_ratio);

    // Define the methods for iter_elastic, iter_inelastic,
    std::vector<std::pair<Target*, Process*>> iter_elastic();
    std::vector<std::pair<Target*, Process*>> iter_inelastic();
    std::vector<std::pair<Target*, Process*>> iter_growth();
    std::vector<std::pair<Target*, Process*>> iter_all();
    std::vector<std::pair<Target*, Process*>> iter_momentum();

    void init();
    void reinitTgas(); // Reinitialization of only Tgas changes
    void reinitEN();   // Reinitialization of only E/N changes

    Eigen::VectorXd maxwell(double kT);
    Eigen::VectorXd iterate(const Eigen::VectorXd& f0, double delta = 1e14);
    Eigen::VectorXd converge(Eigen::VectorXd f0, int maxn = 100, double rtol = 1e-5, double delta0 = 1e14, double m = 4.0, bool full = false);

    double rate(const Eigen::VectorXd& F0, Process* k, bool weighted = false);
    double rate(const Eigen::VectorXd& F0, const std::string& k_str, bool weighted = false);
    Process* search(const std::string& signature, const std::string& product);
    Process* search(const std::string& signature);

    double calcu_mobility(const Eigen::VectorXd& F0);     // mobility*N, [1/m/V/s]
    double mobility() const { return m_mobility; }; // mobility * N
    double calcu_diffusivity(const Eigen::VectorXd& F0);    // diffusion coefficient * N, [1/m/s]
    double diffusivity() const {return m_diffusivity; };
    double mean_energy(const Eigen::VectorXd& F0);  // in eV
    double get_Te() {return m_Te; }; // get electron temperature in Kelvin
    double elec_power() const { return m_elec_power; };  // power/N [eV m^3/s]

private:
    std::unique_ptr<Grid> m_grid;

    //std::unordered_map<std::string, Target*> target;
    std::unordered_map<std::string, std::unique_ptr<Target>> m_targetsMap; // smart pointer

    // Set Linear solver with preconditioner
    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteCholesky<double>> solver;
    //Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

    double EN;
    double kT; // gas temperature in eV
    Eigen::VectorXd density;
    Eigen::VectorXd benergy;
    Eigen::VectorXd cenergy;
    Eigen::VectorXd denergy;
    Eigen::VectorXd denergy32;

    Eigen::VectorXd logSlopeF0;

    int n;

    // Variables used in the init function
    Eigen::VectorXd sigma_eps;
    Eigen::VectorXd sigma_m;
    Eigen::VectorXd sigma_mc;
    Eigen::VectorXd W;
    Eigen::VectorXd DA;
    Eigen::VectorXd DB;

    double _norm(const Eigen::VectorXd& f);
    Eigen::VectorXd _normalized(const Eigen::VectorXd& f);
    Eigen::VectorXd _g(const Eigen::VectorXd& F0);
    Eigen::SparseMatrix<double> _PQ(const Eigen::VectorXd& F0, std::vector<std::pair<Target*, Process*>>& reactions);
    Eigen::SparseMatrix<double> _PQ(const Eigen::VectorXd& F0); // Overloaded version
    Eigen::SparseMatrix<double> _scharf_gummel(const Eigen::VectorXd& sigma_tilde, Eigen::VectorXd G);
    std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> _linsystem(const Eigen::VectorXd& F);

    std::unordered_map<std::string, Process*> m_search_cache;

    // Some important output parameters
    double m_Te;  // electron temperature in K
    double m_mobility; // mobility*N, [1/m/V/s]
    double m_diffusivity; // diffusion coefficient *N, [1/m/s]
    double m_elec_power; // Power/N (eV m^3/s) absorbed by electrons from E-field

    // Cached results for fast calculation of mobility and diffusivity
    Eigen::VectorXd m_DF0;
    std::vector<std::pair<Target*, Process*>> m_growReactions;



};

}
#endif // SOLVER_H

