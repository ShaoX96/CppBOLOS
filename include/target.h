#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include <map>
#include <memory> // for shared_ptr or unique_ptr if we use them for Process objects
#include "process.h"

namespace CppBOLOS {

class Target {
    friend class BoltzmannSolver;
public:
    // Constructor
    Target(const std::string& name);

    // Member functions
    void add_process(Process& process);
    void ensure_elastic();

    // Accessor methods
    std::string getName() const;

    // These might return combined lists of different processes, similar to the Python properties
    std::vector<Process>& inelastic();
    std::vector<Process> everything();

    // Overloaded operators for string representations (if necessary)
    friend std::ostream& operator<<(std::ostream& os, const Target& target);

    double density;

private:
    // Member variables
    std::string name;
    double mass_ratio;  // Initialized to some sentinel value (e.g., -1) to represent "undefined"

    // Lists to store different types of processes
    std::vector<Process> elastic;
    std::vector<Process> effective;
    std::vector<Process> attachment;
    std::vector<Process> ionization;
    std::vector<Process> excitation;
    std::vector<Process> weighted_elastic;
    std::vector<Process> combined_inelastic;

    // Map from process type names to the respective lists (can be implemented differently in C++)
    std::map<std::string, std::vector<Process>*> kind;

    // Map from a product to a list of processes that produce it
    // Prefer not to use pointer in case process.product doesn't exist
    std::map<std::string, std::vector<Process>> by_product;

    // Might need other member variables or utility functions
};

}

#endif //TARGET_H
