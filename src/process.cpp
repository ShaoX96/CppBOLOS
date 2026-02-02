#include "process.h"
#include "utilities.h"
#include <iostream>
#include "Logger.h"

namespace CppBOLOS {

const std::map<std::string, int> Process::IN_FACTOR = {
        {"EXCITATION", 1},
        {"IONIZATION", 2},
        {"ATTACHMENT", 0},
        {"ELASTIC", 1},
        {"MOMENTUM", 1},
        {"EFFECTIVE", 1}
};

const std::map<std::string, int> Process::SHIFT_FACTOR = {
        {"EXCITATION", 1},
        {"IONIZATION", 2},
        {"ATTACHMENT", 1},
        {"ELASTIC", 1},
        {"MOMENTUM", 1},
        {"EFFECTIVE", 1}
};

// Constructor definition in the Process class
Process::Process(
            std::string target,
            std::string kind,
            std::vector<std::vector<double>> data,
            std::string comment,
            double mass_ratio,
            std::string product,
            double threshold,
            double weight_ratio
        )
        : target_name(target), kind(kind), data(data), comment(comment),
          mass_ratio(mass_ratio), product(product), threshold(threshold),
          weight_ratio(weight_ratio), isnull(false), cached_grid(), target(nullptr)
{
    // Copying the 2D data vector into x and y
    for (auto& row : data) {
        this->x.push_back(row[0]);
        this->y.push_back(row[1]);
    }

    // Call the function padinterp and assign the returned function to interp
    if (data.empty()) { LOG_WARNING("CppBOLOS::Process constructor. Empty data input for interp");}
    interp = padinterp(data);

    // Access IN_FACTOR and SHIFT_FACTOR
    auto in_it = IN_FACTOR.find(kind);
    in_factor = (in_it != IN_FACTOR.end()) ? in_it->second : 0;

    auto shift_it = SHIFT_FACTOR.find(kind);
    shift_factor = (shift_it != SHIFT_FACTOR.end()) ? shift_it->second : 0;

    // Check for negative energy
    if (!x.empty() && *std::min_element(x.begin(), x.end()) < 0) {
        throw std::invalid_argument("CppBOLOS::Process: Negative energy in the cross section.");
    }

    // Check for negative cross section
    if (!y.empty() && *std::min_element(y.begin(), y.end()) < 0) {
        throw std::invalid_argument("CppBOLOS::Process: Negative cross section.");
    }
}

Eigen::VectorXd Process::scatterings(const Eigen::VectorXd& g, const Eigen::VectorXd& eps) {
    Eigen::VectorXd r(this->eps.size());

    for (size_t i = 0; i < this->eps.size(); ++i) {
        int idx = this->j[i];
        double gj = g[idx];
        double epsj = eps[idx];
        r[i] = int_linexp0(this->eps[i][0], this->eps[i][1],
                           this->sigma[i][0], this->sigma[i][1],
                           gj, epsj);
    }

    return r;
}

void Process::set_grid_cache(const Grid& grid) {

    if (&grid == cached_grid) {
        // Grid hasn't changed, no need to recompute
        return;
    }

    cached_grid = &grid;

    Eigen::VectorXd eps1 = shift_factor * grid.get_b().array() + threshold;
    eps1 = eps1.array().max((grid.get_b()[0] + 1e-9));
    eps1 = eps1.array().min((grid.get_b().tail<1>()[0] - 1e-9));

    std::vector<double> nodes;
    for (const auto& e : eps1) {
        if (std::find(nodes.begin(), nodes.end(), e) == nodes.end()) {
            nodes.push_back(e);
        }
    }

    for (const auto& b : grid.get_b()) {
        if (b >= eps1[0] && b <= eps1[eps1.size()-1] && std::find(nodes.begin(), nodes.end(), b) == nodes.end()) {
            nodes.push_back(b);
        }
    }

    for (const auto& e : this->x) {
        if (e >= eps1[0] && e <= eps1[eps1.size()-1] && std::find(nodes.begin(), nodes.end(), e) == nodes.end()) {
            nodes.push_back(e);
        }
    }

    std::sort(nodes.begin(), nodes.end());

    std::vector<double> sigma0(nodes.size());
    for (size_t k = 0; k < nodes.size(); ++k) {
        sigma0[k] = this->interp(nodes[k]);
    }

    this->j.resize(nodes.size() - 1);
    this->i.resize(nodes.size() - 1);

    for (int idx = 1; idx < nodes.size(); ++idx) {
        auto it0 = std::lower_bound(grid.b.data(), grid.b.data() + grid.b.size(), nodes[idx]);
        this->j[idx - 1] = std::distance(grid.b.data(), it0) - 1;
        auto it1 = std::lower_bound(eps1.data(), eps1.data() + eps1.size(), nodes[idx]);
        this->i[idx - 1] = std::distance(eps1.data(), it1) - 1;
    }

    this->sigma.resize(nodes.size() - 1, std::vector<double>(2));
    this->eps.resize(nodes.size() - 1, std::vector<double>(2));
    for (size_t idx = 0; idx < nodes.size() - 1; ++idx) {
        this->sigma[idx][0] = sigma0[idx];
        this->sigma[idx][1] = sigma0[idx + 1];
        this->eps[idx][0] = nodes[idx];
        this->eps[idx][1] = nodes[idx + 1];
    }
}

std::ostream& operator<<(std::ostream& os, const Process& process) {
    os << "{" << process.kind << ": " << process.target_name;
    if (!process.product.empty()) {
        os << " -> " << process.product;
    }
    os << "}";
    return os;
}

std::function<double(double)> Process::padinterp(const std::vector<std::vector<double>>& data) const {
    if (data.empty()) {
        // Return an interpolation function that always returns 0
        return [](double) { return 0.0; };
    }

    // Create x and y vectors for interpolation
    std::vector<double> x_interp;
    std::vector<double> y_interp;

    // Prepend zero if needed
    if (data[0][0] > 0) {
        x_interp.push_back(0.0);
        y_interp.push_back(data[0][1]);
    }

    // Copy data
    for (const auto& row : data) {
        x_interp.push_back(row[0]);
        y_interp.push_back(row[1]);
    }

    // Append large value for extrapolation
    x_interp.push_back(1e8);
    y_interp.push_back(data.back()[1]);

    Eigen::VectorXd xXd = Eigen::Map<const Eigen::VectorXd>(x_interp.data(), x_interp.size());
    Eigen::VectorXd yXd = Eigen::Map<const Eigen::VectorXd>(y_interp.data(), y_interp.size());

    // Capture by value to avoid dangling references
    return [xXd, yXd](double xi) {
        return ScipyUtils::interp1d(xXd, yXd, xi);
    };
}


// Integration of u(x) * exp(g * (x0 - x)) * x in [a, b]
// assuming that u is linear with u({a, b}) = {u0, u1}.
inline double Process::int_linexp0(double a, double b, double u0, double u1, double g, double x0) const{
    double result;

    double c0 = (a * u1 - b * u0) / (a - b);
    double c1 = (u0 - u1) / (a - b);

    // When g is zero, we are integrating u(x)*x from a to b
    if (std::abs(g) < 1e-16) {
        // LOG_DEBUG("Process::int_linexp0: |g| < 1e-16.");
        result = c0 * (b*b - a*a) / 2.0 + c1 * (b*b*b - a*a*a) / 3.0;
        return result;
    }

    // Compute the various quantities needed for the calculation
    double expm1a = std::expm1(g * (-a + x0));
    double expm1b = std::expm1(g * (-b + x0));

    double ag = a * g;
    double bg = b * g;

    double ag1 = ag + 1.0;
    double bg1 = bg + 1.0;

    double g2 = g * g;
    double g3 = g2 * g;

    double A1 = (  expm1a * ag1 + ag
                   - expm1b * bg1 - bg) / g2;

    double A2 = (expm1a * (2.0 * ag1 + ag * ag) + ag * (ag + 2.0) -
                 expm1b * (2.0 * bg1 + bg * bg) - bg * (bg + 2.0)) / g3;

    // Compute the result and store it in the memory pointed to by r
    result = c0 * A1 + c1 * A2;

    // When either F0 or F1 is 0 we return 0
    return (std::isnan(result) ? 0.0 : result);
}

// Constructor definition in the NullProcess class
NullProcess::NullProcess(std::string target, std::string kind)
        : Process(target, kind, std::vector<std::vector<double>>{}, "", 0.0, "", 0.0, 0.0) {
    this->isnull = true;
}

// Overloaded output operator definition for NullProcess
std::ostream& operator<<(std::ostream& os, const NullProcess& nullProcess) {
    os << "{NULL}";
    return os;
}

}