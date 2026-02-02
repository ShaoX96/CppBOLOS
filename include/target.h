#ifndef TARGET_H
#define TARGET_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include "process.h"

namespace CppBOLOS {

class Target {
    friend class BoltzmannSolver;
public:
    // Constructor
    Target(const std::string& name);

    // Member functions
    void add_process(std::unique_ptr<Process> process);
    void ensure_elastic();

    // Return combined lists of different processes
    std::vector<Process*> inelastic();
    std::vector<Process*> everything();

    // Overloaded operators for string representations
    friend std::ostream& operator<<(std::ostream& os, const Target& target);

    std::string name;
    double density = 0.0;
    double mass_ratio = -1;  // Initialized to some sentinel value (e.g., -1) to represent "undefined"

private:
    // Lists to store different types of processes
    std::vector<std::unique_ptr<Process>> elastic;
    std::vector<std::unique_ptr<Process>> effective;
    std::vector<std::unique_ptr<Process>> attachment;
    std::vector<std::unique_ptr<Process>> ionization;
    std::vector<std::unique_ptr<Process>> excitation;
    std::vector<std::unique_ptr<Process>> weighted_elastic;

    std::vector<Process*> combined_inelastic;

    // Map from process type names to the respective lists
    std::map<std::string, std::vector<std::unique_ptr<Process>>*> kind;

    // Map from a product to a list of processes that produce it
    std::map<std::string, std::vector<Process*>> by_product;

};

}

#endif //TARGET_H
