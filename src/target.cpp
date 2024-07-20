#include "target.h"
#include <iostream>
#include <sstream>
#include "Logger.h"

namespace CppBOLOS {

Target::Target(const std::string& name)
        : name(name),
          mass_ratio(-1.0),   // We use -1.0 as a sentinel value to represent "undefined"
          density(0.0)
{
    // Initialize the kind map
    kind["ELASTIC"] = &elastic;
    kind["EFFECTIVE"] = &effective;
    kind["MOMENTUM"] = &effective;
    kind["ATTACHMENT"] = &attachment;
    kind["IONIZATION"] = &ionization;
    kind["EXCITATION"] = &excitation;
    kind["WEIGHTED_ELASTIC"] = &weighted_elastic;

    LOG_INFO("Target " + name + " created.")

    // Note: The by_product map in C++ will be empty by default and will be populated as processes are added.
}

void Target::add_process(Process& process) {
    // Retrieve the appropriate kind vector and add the process
    kind[process.kind]->push_back(process);

    // Check and update the mass_ratio if necessary
    if (process.mass_ratio >= 0.0) {  // Assuming -1.0 indicates an undefined mass ratio
        LOG_DEBUG("Mass ratio "  << process.mass_ratio << " for " << this->name);
        if (this->mass_ratio >= 0.0 && this->mass_ratio != process.mass_ratio) {
            throw std::runtime_error("CppBOLOS::add_processError: More than one mass ratio for target '" + this->name + "'");
        }
        this->mass_ratio = process.mass_ratio;
    }

    // Set the process's target
    process.target = this;
    // Update the by_product map
    by_product[process.product].push_back(process);

    LOG_DEBUG("Process " << process << " added to target " << this->name);
}

void Target::ensure_elastic() {
    // Check if both elastic and effective cross-sections are defined
    if (!elastic.empty() && !effective.empty()) {
        throw std::runtime_error("CppBOLOS::ensure_elasticError: In target '" + this->name
        + "': EFFECTIVE and ELASTIC cross-sections are incompatible.");
    }

    // If elastic cross-sections exist, we're done
    if (!elastic.empty()) {
        LOG_DEBUG("Elastic cross-sections already exist.");
        return;
    }

    // If there's more than one effective cross-section, it's an error
    if (effective.size() > 1) {
        throw std::runtime_error("CppBOLOS::ensure_elasticError: In target '" + this->name
        + "': Can't handle more than 1 EFFECTIVE for a given target.");
    }

    // If no effective cross-section exists, log a warning
    if (effective.empty()) {
        LOG_WARNING("Target " + this->name + " has no ELASTIC or EFFECTIVE cross sections")
        return;
    }

    // Extract data from the effective process
    std::vector<std::vector<double>> newdata = effective[0].data;; // Copy of effective data

    // For each inelastic process, subtract from newdata
    for (const Process& p : inelastic()) { // Assuming inelastic() returns combined list of inelastic processes
        for (size_t i = 0; i < newdata.size(); i++) {
            double energy = newdata[i][0];
            newdata[i][1] -= p.interp(energy); // Subtract the interpolated value
        }
    }

    // Check for negative values and set to zero
    bool hasNegative = false;
    for (auto& row : newdata) {
        if (row[1] < 0) {
            row[1] = 0.0;
            hasNegative = true;
        }
    }
    if (hasNegative) {
        LOG_WARNING("After subtracting INELASTIC from EFFECTIVE, target " + this->name +
        " has negative cross-section. Setting as max(0, ...)")
    }
    // Create a new elastic process from the computed data
    Process newElastic(this->name, "ELASTIC", newdata, "Calculated from EFFECTIVE cross sections", effective[0].mass_ratio);

    // Add the new elastic process to this target
    add_process(newElastic);

    LOG_DEBUG("EFFECTIVE -> ELASTIC for target " + this->name)

    // Clear the effective processes
    effective.clear();
}

std::vector<Process>& Target::inelastic(){
    if (combined_inelastic.empty()) {
        combined_inelastic = this->attachment;
        combined_inelastic.insert(combined_inelastic.end(), this->ionization.begin(), this->ionization.end());
        combined_inelastic.insert(combined_inelastic.end(), this->excitation.begin(), this->excitation.end());
    }
    return combined_inelastic; // returns reference to a std::vector<Process> Avoid creating temporary objects.
}

std::vector<Process> Target::everything() {
    std::vector<Process> all_processes = this->elastic;
     std::vector<Process>& inelastic_processes = this->inelastic();
    all_processes.insert(all_processes.end(), inelastic_processes.begin(), inelastic_processes.end());
    return all_processes;
}

std::ostream& operator<<(std::ostream& os, const Target& target) {
    os << target.name;
    return os;
}

}