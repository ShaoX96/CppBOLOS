#include <iostream>
#include <vector>
#include "Eigen/Dense"
#include "Eigen/Sparse"

namespace ScipyUtils {

    // For example, a function equivalent to scipy's interp1d
    double interp1d(const Eigen::VectorXd& x, const Eigen::VectorXd& y, double xi);

    double simpsons_rule(const Eigen::VectorXd& f_values, const Eigen::VectorXd& grid_points);

    // Any other scipy functions you need equivalents for...

}

// Assuming trim is defined something like this:
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (std::string::npos == first) {
        return str;
    }
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// Custom helper function for debugging purpose

void vecPrint(const std::vector<double>& vec, const bool printAll = true, size_t threshold = 5);

void vecPrint(const Eigen::VectorXd& vec, const bool printAll = true, size_t threshold = 5);

void checkMatrixDiag(const Eigen::SparseMatrix<double>& Q, const size_t n = 5);

void checkMatrix(const Eigen::SparseMatrix<double>& Q, const size_t n = 5);

std::vector<double> vecConvert(const Eigen::VectorXd& v);

Eigen::VectorXd vecConvert(const std::vector<double>& v);