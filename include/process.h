#ifndef PROCESS_H
#define PROCESS_H

#include <vector>
#include <map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <functional>
#include "grid.h"
#include "Eigen/Dense"

namespace CppBOLOS {

class Target; // Forward declaration for Target class

class Process {
    friend class Target;
    friend class BoltzmannSolver;
public:
    // The factor of in-scattering.
    static const std::map<std::string, int> IN_FACTOR;

    // The shift factor for inelastic collisions.
    static const std::map<std::string, int> SHIFT_FACTOR;

    Process(
            std::string target = "",
            std::string kind = "",
            std::vector<std::vector<double>> data = {},
            std::string comment = "",
            double mass_ratio = -1.0,
            std::string product = "",
            double threshold = 0.0,
            double weight_ratio = -1.0
    );

    // Delete copy constructor and copy assignment operator
    Process(const Process&) = delete;
    Process& operator=(const Process&) = delete;

    Eigen::VectorXd scatterings(const Eigen::VectorXd& g, const Eigen::VectorXd& eps);

    void set_grid_cache(const Grid& grid);

    friend std::ostream& operator<<(std::ostream& os, const Process& process);

    // Helper functions.
    std::function<double(double)> padinterp(const std::vector<std::vector<double>>& data) const;

    double int_linexp0(double a, double b, double u0, double u1, double g, double x0) const;

    const double get_threshold() const {return threshold;}

private:
    Target* target;

    std::string target_name;
    std::string kind;
    std::string comment;
    std::string product;

    double mass_ratio;
    double threshold;
    double weight_ratio;

    int in_factor;
    int shift_factor;

    std::vector<std::vector<double>> data;
    std::vector<double> x;
    std::vector<double> y;

    const Grid* cached_grid = nullptr;  // Initialize to nullptr as default

    std::vector<int> j;
    std::vector<int> i;
    std::vector<std::vector<double>> sigma;
    std::vector<std::vector<double>> eps;

    std::function<double(double)> interp;

protected:
    bool isnull;

};

class NullProcess : public Process {
public:
    // Constructor
    NullProcess(std::string target, std::string kind);

    // Overloaded output operator declaration
    friend std::ostream& operator<<(std::ostream& os, const NullProcess& nullProcess);
};

}

#endif // PROCESS_H
