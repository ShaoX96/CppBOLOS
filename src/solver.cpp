#include "solver.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include "utilities.h"
#include <stdexcept>
#include "Logger.h"

namespace CppBOLOS {

LogLevel currentLogLevel = LOG_DEBUG;

void BoltzmannSolver::set_grid(std::unique_ptr<Grid> grid) {

    m_grid = std::move(grid);

    this->benergy = m_grid->b;

    this->cenergy = m_grid->c;

    this->denergy = m_grid->d;

    // This is useful when integrating the growth term.
    this->denergy32.resize(this->benergy.size() - 1);

    for (size_t i = 0; i < this->benergy.size() - 1; ++i) {
        this->denergy32[i] = std::pow(this->benergy[i + 1], 1.5) - std::pow(this->benergy[i], 1.5);
    }

    this->n = m_grid->n;

    m_DF0.resize(this->n + 1);
    // Initilize m_DF0, first and last element set to 0.0
    m_DF0(0) = 0.0;
    m_DF0(this->n) = 0.0;
}

void BoltzmannSolver::set_grid(const std::string& type, double start, double end, int N_cells){
    if (type == "LinearGrid") {
        std::unique_ptr<LinearGrid> linGrid(new LinearGrid(start, end, N_cells));
        set_grid(std::move(linGrid));
    } else if (type == "LogGrid") {
        std::unique_ptr<LogGrid> logGrid(new LogGrid(start, end, N_cells));
        set_grid(std::move(logGrid));
    } else if (type == "QuadraticGrid"){
        std::unique_ptr<QuadraticGrid> quadGrid(new QuadraticGrid(start, end, N_cells));
        set_grid(std::move(quadGrid));
    } else if (type == "GeometricGrid"){
        std::unique_ptr<GeometricGrid> geoGrid(new GeometricGrid(start, end, N_cells));
        set_grid(std::move(geoGrid));
    } // TODO: AutomaticGrid
    else {
        throw std::runtime_error("CppBOLOS: Unknown grid type!");
    }
}

void BoltzmannSolver::set_density(const std::string& species, const double& density) {
    if (m_targetsMap.find(species) == m_targetsMap.end()) {
        LOG_WARNING("Species not found in targets map: " + species +", adding new...")
        m_targetsMap[species] = std::make_unique<Target>(species);
    }

    m_targetsMap[species]->density = density;
}

void BoltzmannSolver::set_density(const std::string& ss) {
    std::istringstream stream(ss); // Create a stream from the input string
    std::string token;             // Hold each species:density pair

    while (std::getline(stream, token, ',')) {  // Split the string by comma
        std::string species;
        double dens;
        std::istringstream pairStream(token);   // Stream to parse each species:density pair
        std::getline(pairStream, species, ':'); // Split the pair by colon
        pairStream >> dens;                     // Read the remaining part as density

        // Trim leading and trailing white spaces from species
        species.erase(species.begin(), std::find_if(species.begin(), species.end(), [](int ch) {
            return !std::isspace(ch);
        }));
        species.erase(std::find_if(species.rbegin(), species.rend(), [](int ch) {
            return !std::isspace(ch);
        }).base(), species.end());

        // Now call the set_density function with species and density
        set_density(species, dens);
    }
}

void BoltzmannSolver::set_density(const std::map<std::string, double>& xMap) {
    for (const auto& [sp, dens] : xMap) {
        set_density(sp, dens);
    }

    // Prior to C++17:
    // for (const auto& pair : xMap) {
    //     set_density(pair.first, pair.second);
    // }
}

void BoltzmannSolver::load_collisions(const std::vector<Collision>& collisions) {
    // Loop over the list of collisions and add each one
    for (const auto& collision : collisions) {
        const std::string& target_name = collision.target;

        // Try to get the target from the map
        if (m_targetsMap.find(target_name) == m_targetsMap.end()) {
            // If the target doesn't exist, create a new one and add it to the map
            m_targetsMap[target_name] = std::make_unique<Target>(target_name);
        }

        Target* target = m_targetsMap[target_name].get();

        // Create a new Process
        auto process = std::make_unique<Process>(
                target_name,
                collision.kind,
                collision.data,
                "",
                collision.mass_ratio,
                collision.product,
                collision.threshold,
                collision.weight_ratio
        );

        // Add the process to the target
        target->add_process(std::move(process));
    }

    // Ensure all targets have their elastic cross-sections
    for (auto& [key, targetPtr] : m_targetsMap) {
        targetPtr->ensure_elastic();
    }

    // Extra initialization for growthProcesses
    m_growReactions = iter_growth();
    for (auto& [target, process] : m_growReactions) {
        process->set_grid_cache(*m_grid);
    }
}

std::vector<std::pair<Target*, Process*>> BoltzmannSolver::iter_elastic() {
    std::vector<std::pair<Target*, Process*>> results;

    for (const auto& [key, targetPtr] : m_targetsMap) {
        Target* target = targetPtr.get();
        if (target->density > 0) {
            for (const auto& processPtr : target->elastic) {
                results.emplace_back(target, processPtr.get());
            }
        }
    }
    return results;
}

std::vector<std::pair<Target*, Process*>> BoltzmannSolver::iter_inelastic() {
    std::vector<std::pair<Target*, Process*>> results;

    for (auto& [key, targetPtr] : m_targetsMap) {
        Target* target = targetPtr.get();
        if (target->density > 0) {
            for (Process* process : target->inelastic()) {
                results.emplace_back(target, process);
            }
        }
    }

    return results;
}

std::vector<std::pair<Target*, Process*>> BoltzmannSolver::iter_growth() {
    std::vector<std::pair<Target*, Process*>> results;

    for (auto& [key, targetPtr] : m_targetsMap) {
        Target* target = targetPtr.get();
        if (target->density > 0) {
            for (const auto& processPtr : target->ionization) {
                results.emplace_back(target, processPtr.get());
            }
            for (const auto& processPtr : target->attachment) {
                results.emplace_back(target, processPtr.get());
            }
        }
    }

    return results;
}

std::vector<std::pair<Target*, Process*>> BoltzmannSolver::iter_all() {
    std::vector<std::pair<Target*, Process*>> results;

    // results.reserve(this->target.size() * 2);  // Assuming on average 2 processes per target, adjust as necessary.

    // Directly insert elastic processes into results
    for (const auto& procPair : iter_elastic()) {
        results.push_back(procPair);
    }

    // Directly insert inelastic processes into results
    for (const auto& procPair : iter_inelastic()) {
        results.push_back(procPair);
    }

    return results;
}

std::vector<std::pair<Target*, Process*>> BoltzmannSolver::iter_momentum() {
    return iter_all();
}

void BoltzmannSolver::init() {
    // Initialize sigma_eps and sigma_m with zeros
    this->sigma_eps.setZero(this->benergy.size()); // Equivalent to np.zeros_like
    this->sigma_m.setZero(this->benergy.size());
    sigma_mc.setZero(this->cenergy.size());

    // Iterating over elastic processes
    std::vector<std::pair<Target*, Process*>> elasticProcesses = iter_elastic();
    for (auto& [target, process] : elasticProcesses) {
        Eigen::VectorXd s(this->benergy.size());
        for (size_t i = 0; i < s.size(); i++) {
            s[i] = target->density * process->interp(this->benergy[i]);
        }
        this->sigma_eps += 2 * target->mass_ratio * s;
        this->sigma_m += s;
        process->set_grid_cache(*m_grid);
    }

    // Iterating over inelastic processes
    std::vector<std::pair<Target*, Process*>> inelasticProcesses = iter_inelastic();
    for (auto& [target, process] : inelasticProcesses) {
        if (!process) {
            LOG_DEBUG("Invalid process pointer encountered!");
            continue;
        }

        Eigen::VectorXd s(this->benergy.size());

        for (size_t i = 0; i < s.size(); i++) {
            s[i] = target->density * process->interp(this->benergy[i]);
        }
        this->sigma_m += s;

        process->set_grid_cache(*m_grid);
    }

    for (const auto& [target, process] : this->iter_momentum()) {
        Eigen::VectorXd s(this->cenergy.size());
        for (int i = 0; i < this->cenergy.size(); ++i) {
            s(i) = target->density * process->interp(this->cenergy(i));
        }
        sigma_mc += s;
    }

    this->W = -GAMMA * this->benergy.array().square() * this->sigma_eps.array();
    this->DA = (GAMMA / 3.0) * std::pow(this->EN, 2) * this->benergy.array();
    this->DB = GAMMA * this->kT * this->benergy.array().square() * this->sigma_eps.array();

    LOG_INFO("Solver successfully initialized/updated.");
}

void BoltzmannSolver::reinitEN() {
    this->DA = (GAMMA / 3.0) * std::pow(this->EN, 2) * this->benergy.array();
    LOG_INFO("Reinitialization for E/N update only");
}

void BoltzmannSolver::reinitTgas() {
    this->DB = GAMMA * this->kT * this->benergy.array().square() * this->sigma_eps.array();
    LOG_INFO("Reinitialization for temperature update only")
}

Eigen::VectorXd BoltzmannSolver::maxwell(double kT) {
    // Calculates a Maxwell-Boltzmann distribution function.

    Eigen::VectorXd maxwell_dist(this->cenergy.size());
    std::transform(this->cenergy.begin(), this->cenergy.end(), maxwell_dist.begin(),
                   [kT](double &c) { return 2 * std::sqrt(1 / M_PI) * std::pow(kT, -1.5) * std::exp(-c / kT); });
    return maxwell_dist;
}

double BoltzmannSolver::_norm(const Eigen::VectorXd& f) {
    // Calculate the norm of the function using Simpson's rule for integration
    Eigen::VectorXd integrand(f.size());
    integrand = f.array() * this->cenergy.array().sqrt();
    double result = ScipyUtils::simpsons_rule(integrand, this->cenergy);
    return result;
}

Eigen::VectorXd BoltzmannSolver::_normalized(const Eigen::VectorXd& f) {
    // Normalize the function
    double N = _norm(f);
    return f / N;
}

// Calculate (local) logarithmic slope. Eq.51 in Ref. Hagelaar 2005
Eigen::VectorXd BoltzmannSolver::_g(const Eigen::VectorXd& F0) {
    // Create Fp (with two additional elements)
    int size = F0.size();
    Eigen::VectorXd Fp(size + 2);
    Fp << F0[0], F0, F0[size - 1];

    // Create cenergyp (with two additional elements)
    Eigen::VectorXd cenergyp(size + 2);
    cenergyp << this->cenergy(0), this->cenergy, this->cenergy(size - 1);

    // Calculate g
    Eigen::VectorXd g(size);
    for (int i = 0; i < size; ++i) {
        g[i] = std::log(Fp[i + 2] / Fp[i]) / (cenergyp[i + 2] - cenergyp[i]);
    }

    return g;
}

Eigen::SparseMatrix<double> BoltzmannSolver::_PQ(const Eigen::VectorXd& F0, std::vector<std::pair<Target*, Process*>>& reactions) {

    Eigen::SparseMatrix<double> PQ(this->n, this->n);

    Eigen::VectorXd g = _g(F0);

    // If no specific reactions are provided, use the iter_inelastic() method
    if (reactions.empty()) {
        reactions = this->iter_inelastic();
    }

    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(reactions.size() * 2 * this->n);

    for (const std::pair<Target*, Process*>& tk : reactions) {
        Process* k = tk.second;
        Target* t = tk.first;

        Eigen::VectorXd r = t->density * GAMMA * k->scatterings(g, this->cenergy);

        double in_factor = k->in_factor;  // Assumes a property or method to get in_factor in C++

        // Assuming i and j vectors are of the same size
        for (size_t idx = 0; idx < k->i.size(); ++idx) {
            tripletList.emplace_back(k->i[idx], k->j[idx], in_factor * r[idx]);
            tripletList.emplace_back(k->j[idx], k->j[idx], -r[idx]);
        }
    }

    PQ.setFromTriplets(tripletList.begin(), tripletList.end());

    return PQ;
}

// Overloaded version of _PQ() which only takes F0
Eigen::SparseMatrix<double> BoltzmannSolver::_PQ(const Eigen::VectorXd& F0) {
    std::vector<std::pair<Target*, Process*>> reactions = this->iter_inelastic();
    return _PQ(F0, reactions);
}

Eigen::SparseMatrix<double> BoltzmannSolver::_scharf_gummel(const Eigen::VectorXd& sigma_tilde, Eigen::VectorXd G) {
    // Calculate D
    Eigen::VectorXd D = this->DA.array() / sigma_tilde.array() + this->DB.array(); // Eq.41

    // Compute the z values
    Eigen::VectorXd z(this->cenergy.size() + 1);
    z[0] = std::nan(""); // Set to NaN
    for (int i = 1; i < this->cenergy.size(); i++) {
        // z[i] = this->W[i] * (this->cenergy[i] - this->cenergy[i - 1]) / D[i];
        if (D[i] == 0) {
            z[i] = std::numeric_limits<double>::infinity();
        } else {
            z[i] = this->W[i] * (this->cenergy[i] - this->cenergy[i - 1]) / D[i];
        }
    }
    z[this->cenergy.size()] = std::nan(""); // Set to NaN

    // Compute a0 and a1 terms
    // Eigen::VectorXd a0 = this->W.array() / (1.0 - (-z.array()).exp());
    // Eigen::VectorXd a1 = this->W.array() / (1.0 - z.array().exp());
    Eigen::VectorXd a0 = (z.array() != 0).select(
            this->W.array() / (1.0 - (-z.array()).exp()),
            std::numeric_limits<double>::infinity()
    );

    Eigen::VectorXd a1 = (z.array() != 0).select(
            this->W.array() / (1.0 - z.array().exp()),
            std::numeric_limits<double>::infinity()
    );

    // Construct the sparse matrix A
    Eigen::SparseMatrix<double> A(this->n, this->n);
    std::vector<Eigen::Triplet<double>> triplets;

    // Populate the diagonals using Eigen's array operations
    // NOTE: carefully check the size of a0, a1
    triplets.reserve(3 * this->n); // Estimation to preallocate space
    triplets.emplace_back(0, 0, a0[1] + G[0]);
    triplets.emplace_back(0, 1, a1[1]);
    for (int i = 1; i < this->n - 1; i++) {
        triplets.emplace_back(i, i, a0[i+1]-a1[i] + G[i]); // Main diagonal
        triplets.emplace_back(i, i - 1, -a0[i]);           // Sub-diagonal
        triplets.emplace_back(i, i + 1, a1[i+1]);          // Super-diagonal
    }

    // Handle the boundary conditions directly using Eigen's segment function
    triplets.emplace_back(this->n - 1, this->n - 2, -a0[this->n - 1]);
    triplets.emplace_back(this->n - 1, this->n - 1, -a1[this->n - 1]);

    A.setFromTriplets(triplets.begin(), triplets.end());

    return A;
}

std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> BoltzmannSolver::_linsystem(const Eigen::VectorXd& F0) {

    // Calculate Q
    Eigen::SparseMatrix<double> Q = _PQ(F0);

    // Calculate nu
    double nu = (Q * F0).sum();

    // Compute sigma_tilde with division by zero handling
    // Eigen::VectorXd sigma_tilde = this->sigma_m.array() + nu / this->benergy.array().sqrt() / GAMMA;
    Eigen::VectorXd sigma_tilde = sigma_m.array() + (benergy.array() != 0).select(
            nu / benergy.array().sqrt() / GAMMA,
            std::numeric_limits<double>::infinity()
    );

    // Compute the R (G) term
    Eigen::VectorXd G = 2.0 * this->denergy32 * nu / 3.0;

    // Compute matrix A
    Eigen::SparseMatrix<double> A = _scharf_gummel(sigma_tilde, G);

    return std::make_pair(A, Q);
}

Eigen::VectorXd BoltzmannSolver::iterate(const Eigen::VectorXd& F0, double delta) {
    // Get matrices A and Q from the _linsystem method
    std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>> matrices = _linsystem(F0);
    Eigen::SparseMatrix<double> A = matrices.first;
    Eigen::SparseMatrix<double> Q = matrices.second;

    // checkMatrix(A);
    // checkMatrix(Q);
    // Set up the linear system
    Eigen::SparseMatrix<double> I(F0.size(), F0.size());  // Identity matrix
    I.setIdentity();
    Eigen::SparseMatrix<double> lhs = I + delta * A - delta * Q;

    // if (!lhs.isApprox(lhs.transpose())) {
    //     std::cerr << "Matrix is not symmetric!" << std::endl;
    // }

    // Solve the system using Eigen's solver. The BiCGSTAB solver is a good choice for sparse systems.
    solver.compute(lhs);
    if(solver.info() != Eigen::Success) {
        throw std::runtime_error("BoltzmannSolver::iterate(): BiCGSTAB solver decomposition failed!");
    }

    Eigen::VectorXd f1 = solver.solve(F0);
    if(solver.info() != Eigen::Success) {
        throw std::runtime_error("BoltzmannSolver::iterate(): BiCGSTAB solving failed!");
    }

    // Ensure f1 has non-zero values
    f1 = f1.unaryExpr([](double v) { return std::max(v, 1E-100); });

    // Return the normalized solution
    return _normalized(f1);
}

Eigen::VectorXd BoltzmannSolver::converge(const Eigen::VectorXd& F0, int maxn, double rtol, double delta0, double m, bool full) {
    double err0 = 0.0, err1 = 0.0;
    double delta = delta0;
    Eigen::VectorXd f1;
    Eigen::VectorXd f0 = F0;

    for (int i = 0; i < maxn; ++i) {
        // Richardson Extrapolation for updating delta
        // if (err1 > 0 && err1 < err0) {
        //     double new_delta = delta * std::log(m) / (std::log(err0) - std::log(err1));
        //     delta = delta + 0.5 * (new_delta - delta);
        // } else {
        //     delta = delta0;
        // }
        // LOG_DEBUG("delta = " << delta);

        // Perform Iteration
        try {
            f1 = this->iterate(f0, delta);
        } catch (const std::runtime_error& e) {
            LOG_DEBUG(e.what())
            LOG_DEBUG("Iterate failed with delta = " << delta << ". Reducing delta and retrying.")
            delta = delta / 2.0;
            f0 = F0;
            continue;
        }

        // Compute Error
        err0 = err1;
        err1 = this->_norm((f0 - f1).cwiseAbs());

        // Logging (this is a simple implementation; you might want to use a proper logging library)
        LOG_DEBUG( "After iteration " << (i + 1) << ", err = " << err1 << " (target: " << rtol << ")");

        // Check Convergence
        if (err1 < rtol) {
            LOG_INFO("***** Convergence achieved after " << (i + 1) << " iterations. err = " << err1 << " *****");
            nIter = i+1;

            if (full) {
                // TODO: Return additional convergence information if needed
                Eigen::VectorXd f1Full(f1.size() + 2);
                f1Full << f1, i + 1, err1;
                return f1Full;
            }

            // ------------ Populate some parameters ------------
            m_Te = 2.0/3.0 * EV2K * mean_energy(f1);

            logSlopeF0 = _g(f1); // used for rate calculation

            for (int i = 1; i < f1.size(); ++i) {
                double diff_F0 = f1(i) - f1(i - 1);
                double diff_cenergy = cenergy(i) - cenergy(i - 1);
                m_DF0(i) = diff_F0 / diff_cenergy;
            }

            // Calculate Q
            Eigen::SparseMatrix<double> Q = this->_PQ(f1, m_growReactions);

            // Calculate nu
            double nu = (Q * f1).sum() / GAMMA;
            // Calculate sigma_tilde
            // Eigen::VectorXd sigma_tilde = this->sigma_m.array() + nu / this->benergy.array().sqrt();
            Eigen::VectorXd sigma_tilde = sigma_m.array() + (benergy.array() != 0).select(
                    nu / benergy.array().sqrt(),
                    std::numeric_limits<double>::infinity()
            );

            // Calculate y
            Eigen::VectorXd y = m_DF0.array() * this->benergy.array() / sigma_tilde.array();
            y[0] = 0;

            // Integrate y with respect to benergy to find mobility
            m_mobility = -(GAMMA / 3) * ScipyUtils::simpsons_rule(y, this->benergy);  // Assumes integrate() is a utility function

            // Calculate diffusivity
            sigma_tilde = sigma_mc.array() + (cenergy.array() != 0).select(
                    nu / cenergy.array().sqrt(),
                    std::numeric_limits<double>::infinity()
            );

            y = f1.array() * this->cenergy.array() / sigma_tilde.array();

            // Integrate y with respect to cenergy to find diffusion coefficient
            m_diffusivity = (GAMMA / 3) * ScipyUtils::simpsons_rule(y, this->cenergy);  // Assumes integrate() is a utility function

            m_elec_power = m_mobility * EN * EN;

            return f1;
        }

        f0 = f1;
    }

    // Convergence was not achieved
    throw ConvergenceError("Convergence failed after " + std::to_string(maxn) + " iterations!");
}

// Search for a process or a number of processes within the solver
// If `first` is true, return only the first process matching the search
// Otherwise, return a vector of matches (possibly empty)
Process* BoltzmannSolver::search(const std::string& signature, const std::string& product) {
    // Use signature + product as the cache key
    std::string cache_key = signature + "->" + product;

    // Check if it's already in the cache
    auto it = m_search_cache.find(cache_key);
    if (it != m_search_cache.end()) {
        return it->second; // Return the cached result
    }

    // If product is specified
    if (!product.empty()) {
        auto target_it = m_targetsMap.find(signature);
        if (target_it == m_targetsMap.end()) {
            throw std::runtime_error("CppBOLOS searchError: Target " + signature + " not found");
        }

        // Look for the product in the by_product map
        auto process_it = target_it->second->by_product.find(product);
        if (process_it == target_it->second->by_product.end()) {
            throw std::runtime_error("CppBOLOS searchError: Process for " + signature + " -> " + product + " not found");
        }

        // Return the first process matching the search
        if (! process_it->second.empty()) {
            Process* result = process_it->second[0];
            // Cache the result
            m_search_cache[cache_key] = result;

            return result;
        }
    }
    // Only the signature is specified
    else {
        return search(signature);
    }

    return nullptr;  // Return nullptr if no process is found
}

// Overloaded function to search by signature only
Process* BoltzmannSolver::search(const std::string& signature) {
    // Check if it's already in the cache
    auto it = m_search_cache.find(signature);
    if (it != m_search_cache.end()) {
        return it->second; // Return the cached result
    }

    std::string target_name, product_name;
    size_t arrowPos = signature.find("->");
    if (arrowPos != std::string::npos) {
        target_name = signature.substr(0, arrowPos);
        product_name = signature.substr(arrowPos + 2);

        // Trim whitespaces (if needed)
        target_name = trim(target_name);
        product_name = trim(product_name);
    }

    Process* result = search(target_name, product_name);
    // Cache the result
    m_search_cache[signature] = result;

    return result;
}

double BoltzmannSolver::rate(const Eigen::VectorXd& F0, Process* k, bool weighted) {
    // Set the grid cache for the process k
    k->set_grid_cache(*m_grid);

    // Calculate scatterings
    // TODO: if EEDF is updated from ML technique, recalculate lgoSlopeF0 here
    Eigen::VectorXd r = k->scatterings(logSlopeF0, this->cenergy);

    // Build a sparse matrix for P
    Eigen::SparseMatrix<double> P(this->n, 1);
    std::vector<Eigen::Triplet<double>> tripletList(r.size());

    for (int i = 0; i < r.size(); ++i) {
        tripletList.emplace_back(k->j[i], 0, GAMMA * r[i]);
    }

    P.setFromTriplets(tripletList.begin(), tripletList.end());

    // Convert the sparse matrix to a dense vector and multiply with GAMMA
    Eigen::VectorXd P_dense = Eigen::VectorXd(P);

    // Calculate the rate
    double rate = F0.dot(P_dense);

    // If weighted, multiply by the density of the target
    if (weighted) {
        rate *= k->target->density;  // Assuming target is a pointer to the target object in Process
    }

    return rate;
}

// Overloaded function that takes a string to identify the process
double BoltzmannSolver::rate(const Eigen::VectorXd& F0, const std::string& k_str, bool weighted) {
    // Search for the process using the string identifier (to be implemented)
    Process* k = search(k_str);
    return rate(F0, k, weighted);
}

double BoltzmannSolver::calcu_mobility(const Eigen::VectorXd& F0) {
    // Initialize DF0
    Eigen::VectorXd DF0(F0.size() + 1);
    DF0(0) = 0.0;
    DF0(F0.size()) = 0.0;

    for (int i = 1; i < F0.size(); ++i) {
        double diff_F0 = F0(i) - F0(i - 1);
        double diff_cenergy = cenergy(i) - cenergy(i - 1);
        DF0(i) = diff_F0 / diff_cenergy;
    }

    // Calculate Q
    Eigen::SparseMatrix<double> Q = this->_PQ(F0, m_growReactions);

    // Calculate nu
    double nu = (Q * F0).sum() / GAMMA;

    // Calculate sigma_tilde
    // Eigen::VectorXd sigma_tilde = this->sigma_m.array() + nu / this->benergy.array().sqrt();
    Eigen::VectorXd sigma_tilde = sigma_m.array() + (benergy.array() != 0).select(
            nu / benergy.array().sqrt(),
            std::numeric_limits<double>::infinity()
    );

    // Calculate y
    Eigen::VectorXd y = DF0.array() * this->benergy.array() / sigma_tilde.array();
    y[0] = 0;

    // Integrate y with respect to benergy to find mobility
    double mobility = -(GAMMA / 3) * ScipyUtils::simpsons_rule(y, this->benergy);  // Assumes integrate() is a utility function

    return mobility;
}

// Calculates the diffusion coefficient from a distribution function.
double BoltzmannSolver::calcu_diffusivity(const Eigen::VectorXd& F0) {
    // Calculate Q
    Eigen::SparseMatrix<double> Q = this->_PQ(F0, m_growReactions);

    // Calculate nu
    double nu = (Q * F0).sum() / GAMMA;

    // Calculate sigma_tilde
    // Eigen::VectorXd sigma_tilde = sigma_m.array() + nu / this->cenergy.array().sqrt();
    Eigen::VectorXd sigma_tilde = sigma_mc.array() + (cenergy.array() != 0).select(
            nu / cenergy.array().sqrt(),
            std::numeric_limits<double>::infinity()
    );

    // Calculate y
    Eigen::VectorXd y = F0.array() * this->cenergy.array() / sigma_tilde.array();

    // Integrate y with respect to cenergy to find diffusion coefficient
    double diffusion = (GAMMA / 3) * ScipyUtils::simpsons_rule(y, this->cenergy);  // Assumes integrate() is a utility function

    return diffusion;
}

double BoltzmannSolver::mean_energy(const Eigen::VectorXd& F0) {
    double sum = 0.0;

    // Loop through to calculate the diff(benergy^2.5) for each interval
    for (int i = 0; i < benergy.size() - 1; ++i) {
        double de52 = std::pow(benergy(i + 1), 2.5) - std::pow(benergy(i), 2.5);
        sum += 0.4 * F0(i) * de52;
    }

    return sum;
}

}


