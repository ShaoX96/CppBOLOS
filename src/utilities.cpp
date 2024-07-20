#include "utilities.h"
#include <stack>

namespace ScipyUtils {

double interp1d(const Eigen::VectorXd& x, const Eigen::VectorXd& y, double xi) {
    // This is a basic linear interpolation function.
    // Ensure x and y are sorted.
    auto it = std::lower_bound(x.data(), x.data() + x.size(), xi);
    if (it == x.data() + x.size()) return y[y.size() - 1];
    if (it == x.data()) return y[0];

    int idx = std::distance(x.data(), it);
    double x1 = x[idx - 1];
    double y1 = y[idx - 1];
    double x2 = x[idx];
    double y2 = y[idx];

    return y1 + (xi - x1) * (y2 - y1) / (x2 - x1);
}

double simpsons_rule(const Eigen::VectorXd& f_values, const Eigen::VectorXd& grid_points) {
        double integral = 0.0;
        size_t n = grid_points.size();

        // If the number of points is even, we'll stop one interval earlier and handle the last interval with the Trapezoidal rule.
        size_t stop = (n % 2 == 0) ? n - 2 : n - 1;

        for (size_t i = 0; i < stop; i += 2) {
                double a = grid_points[i];
                double b = grid_points[i + 2];
                double h = (b - a) / 2.0;  // half of the interval size

                integral += (h/3.0) * (f_values[i] + 4*f_values[i + 1] + f_values[i + 2]);
        }

        // If grid_points.size() is even, handle the last interval with the Trapezoidal rule.
        if (n % 2 == 0) {
                double a = grid_points[n - 2];
                double b = grid_points[n - 1];
                integral += (b - a) * 0.5 * (f_values[n - 2] + f_values[n - 1]);
        }

        return integral;
}

    // Implementations of other scipy-equivalent functions...

}

// Custom helper function for debugging purpose
void vecPrint(const std::vector<double>& vec, const bool printAll, size_t threshold) {
    if (vec.empty()) {
        std::cout << "Vector is empty!" << std::endl;
        return;
    }
    std::cout << "-------- Checking Vector --------" << std::endl;

    std::cout << "Size: " << vec.size() << " | First: " << vec[0] << " | Last: " << vec[vec.size() - 1] << std::endl;

    if (printAll) {
        if (vec.size() <= 2*threshold) {
            std::cout << "Vector elements:" << std::endl;
            for (const auto& elem : vec) {
                std::cout << elem << ", ";
            }
            std::cout << std::endl;
        } else {
            std::cout << "First " << threshold << " elements:" << std::endl;
            for (size_t i = 0; i < threshold; ++i) {
                std::cout << vec[i] << ", ";
            }
            std::cout << "..., " << std::endl;
            std::cout << "Last " << threshold << " elements:" << std::endl;
            for (size_t i = vec.size() - threshold; i < vec.size(); ++i) {
                std::cout << vec[i] << ", ";
            }
            std::cout << std::endl;
        }
    }
    std::cout << "---------------------------------" << std::endl;
}

void vecPrint(const Eigen::VectorXd& vec, const bool printAll, size_t threshold){
    std::vector<double> stdVec(vec.data(), vec.data() + vec.size());
    vecPrint(stdVec, printAll, threshold);
}

/*
void checkMatrix(Eigen::SparseMatrix<double> Q, size_t n){
    std::cout << "rows: " << Q.rows() << " | cols: " << Q.cols() << " | first " << n << " data: " << std::endl;
    int printed = 0;
    for(int k = 0; k < Q.outerSize(); ++k) {
        for(Eigen::SparseMatrix<double>::InnerIterator it(Q,k); it; ++it) {
            std::cout << it.value() << std::endl;
            printed++;
            if (printed >= n) { return;}
        }
}
*/

void checkMatrixDiag(const Eigen::SparseMatrix<double>& Q, const size_t n) {
    std::cout << "rows: " << Q.rows() << " | cols: " << Q.cols() << " | first " << n << " diagonals: " << std::endl;

    // Get the minimum between the matrix size and n to avoid out of bounds access
    size_t count = std::min(n, static_cast<size_t>(std::min(Q.rows(), Q.cols())));

    for (size_t i = 0; i < count; ++i) {
        std::cout << Q.coeff(i, i) << " ";
    }
    std::cout << std::endl;
}

void checkMatrix(const Eigen::SparseMatrix<double>& Q, const size_t n) {
    std::cout << "rows: " << Q.rows() << " | cols: " << Q.cols() << std::endl;

    struct MatrixEntry {
        Eigen::Index row;
        Eigen::Index col;
        double value;
    };

    std::vector<MatrixEntry> nonZeroEntries;

    // Storing non-zero entries
    for (int k = 0; k < Q.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(Q, k); it; ++it) {
            nonZeroEntries.push_back({it.row(), it.col(), it.value()});
        }
    }

    // Printing the first 'n' non-zero entries
    for (size_t i = 0; i < std::min(n, nonZeroEntries.size()); ++i) {
        std::cout << "(" << nonZeroEntries[i].row << "," << nonZeroEntries[i].col << "): " << nonZeroEntries[i].value << "  ";
        if ((i + 1) % 5 == 0) {
            std::cout << std::endl;
        }
    }

    std::cout << "...\n";

    // Printing the last 'n' non-zero entries
    for (size_t i = std::max(n, nonZeroEntries.size()) - n; i < nonZeroEntries.size(); ++i) {
        std::cout << "(" << nonZeroEntries[i].row << "," << nonZeroEntries[i].col << "): " << nonZeroEntries[i].value << "   ";
        if ((i + 1) % 5 == 0 || i == nonZeroEntries.size() - 1) {
            std::cout << std::endl;
        }
    }

    std::cout << "\n";
}

// Convert Eigen::VectorXd to std::vector<double>
std::vector<double> vecConvert(const Eigen::VectorXd& v) {
    return {v.data(), v.data() + v.size()};
}

// Convert std::vector<double> to Eigen::VectorXd
Eigen::VectorXd vecConvert(const std::vector<double>& v) {
    return Eigen::VectorXd::Map(v.data(), v.size());
}
