#include "grid.h"
#include "utilities.h"
#include <cmath>
#include "Logger.h"

namespace CppBOLOS {

Grid::Grid(double x0, double x1, int n)
    : x0(x0), x1(x1), n(n), delta(x1 - x0)
{
    // Avoid initializing other member variables here that depend on virtual functions
    _interp = nullptr;
}

void Grid::initialize() {
    fx0 = f(x0);
    fx1 = f(x1);

    Eigen::VectorXd fx = Eigen::VectorXd::LinSpaced(n + 1, fx0, fx1);

    b = fx.unaryExpr([this](double val) { return finv(val); });
    c = 0.5 * (b.head(n) + b.tail(n));

    d = b.segment(1, n) - b.head(n);

    df = fx[1] - fx[0]; // spacing of the mapped x

    //This is useful in some routines that integrate eps**1/2 * f
    d32 = b.segment(1, n).array().pow(1.5) - b.head(n).array().pow(1.5);
}

Eigen::VectorXd Grid::interpolate(const Eigen::VectorXd& f, const Grid& other) {

    if (!_interp) {
        // Prepare x and y values for interpolation
        Eigen::VectorXd x_values(other.c.size() + 2);
        Eigen::VectorXd y_values(f.size() + 2);

        // Populate the x and y vectors
        x_values << other.x0, other.c, other.x1;
        y_values << f[0], f, f[f.size() - 1];

        _interp = [x_values, y_values](double xi) {
            if (xi < x_values(0) || xi > x_values(x_values.size() - 1)) {
                return 0.0;
            }
            return ScipyUtils::interp1d(x_values, y_values, xi);
        };
    }

    Eigen::VectorXd result(this->c.size());
    for (int i = 0; i < this->c.size(); ++i) {
        result(i) = _interp(this->c(i));
    }

    return result;
}

int Grid::cell(double x) const{
    return static_cast<int>((f(x) - fx0) / df);
}

LinearGrid::LinearGrid(double x0, double x1, int n)
    : Grid(x0, x1, n)
{
    // Call the initialization logic after the base class constructor
    initialize();
}

QuadraticGrid::QuadraticGrid(double x0, double x1, int n)
    : Grid(x0, x1, n)
{
    initialize();
}

double QuadraticGrid::f(double x) const{
    return std::sqrt(x - x0);
}

double QuadraticGrid::finv(double w) const{
    return w * w + x0;
}

GeometricGrid::GeometricGrid(double x0, double x1, int n, double r)
    : Grid(x0, x1, n), r(r)
{
    logr = std::log(r);
    rn_minus_1 = std::exp(n * logr) - 1;
    initialize();
}

double GeometricGrid::f(double x) const {
    return (std::log(1 + (x - x0) * rn_minus_1 / delta) / logr);
}

double GeometricGrid::finv(double w) const {
    return (x0 + delta * (std::exp(w * logr) - 1) / rn_minus_1);
}

LogGrid::LogGrid(double x0, double x1, int n, double s)
    : Grid(x0, x1, n), s(s)
{
    initialize();
}

double LogGrid::f(double x) const {
    return std::log(s + (x - x0));
}

double LogGrid::finv(double w) const {
    return std::exp(w) - s + x0;
}

AutomaticGrid::AutomaticGrid(const Grid& grid, const Eigen::VectorXd& f0, double delta)
    : Grid(grid.get_x0(), grid.get_x1(), grid.get_n()), delta(delta)
{
    initialize();

    Eigen::VectorXd cum(n+1);
    for (size_t i = 1; i <= n; i++) {
        cum[i] = cum[i-1] + d32[i-1] * f0[i-1];
    }
    cum /= cum(cum.size() - 1);

    Eigen::VectorXd nnew = Eigen::VectorXd::LinSpaced(n + 1, 0.0, 1.0);
    for (size_t i = 0; i <= n; i++) {
        b[i] = ScipyUtils::interp1d(cum, b, nnew[i]);
    }

    c = 0.5 * (b.segment(1, b.size() - 1) + b.head(b.size() - 1));
    d = b.segment(1, b.size() - 1) - b.head(b.size() - 1);
    _interp = nullptr;

}


/*
std::unique_ptr<Grid> mkgrid(const std::string& kind, ...) { // Again, handling variable arguments is tricky
    if (kind == "linear" || kind == "lin") {
        // Extract relevant arguments and create a LinearGrid object
        return std::make_unique<LinearGrid>(...);
    }
    else if (kind == "quadratic" || kind == "quad") {
        // Extract relevant arguments and create a QuadraticGrid object
        return std::make_unique<QuadraticGrid>(...);
    }
    else if (kind == "geometric" || kind == "geo") {
        // Extract relevant arguments and create a GeometricGrid object
        return std::make_unique<GeometricGrid>(...);
    }
    else if (kind == "logarithmic" || kind == "log") {
        // Extract relevant arguments and create a LogGrid object
        return std::make_unique<LogGrid>(...);
    }
    // ... Add other kinds as needed

    // If the kind does not match any known grid type, return nullptr or throw an exception
    return nullptr;
}
*/

}