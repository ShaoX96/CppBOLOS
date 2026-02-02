#ifndef GRID_H
#define GRID_H

#include <vector>
#include <functional>
#include <string>
#include "Eigen/Dense"

namespace CppBOLOS {

class Grid {
    friend class BoltzmannSolver;
    friend class Process;
public:
    Grid(double x0, double x1, int n);

    //Virtual destructor to ensure proper cleanup for objects of derived types.
    virtual ~Grid() = default;

    Eigen::VectorXd interpolate(const Eigen::VectorXd& f, const Grid& other);
    int cell(double x) const;

    double get_x0() const { return x0; }
    double get_x1() const { return x1; }
    int get_n() const { return n; }
    Eigen::VectorXd get_b() const {return b;}
    Eigen::VectorXd get_c() const {return c;}

protected:
    double x0, x1, delta, fx0, fx1, df;
    int n;
    Eigen::VectorXd b, c, d, d32;
    std::function<double(double)> _interp;

    //Make f and finv pure virtual functions in the Grid class to enforce their implementation in derived classes.
    virtual double f(double x) const = 0;
    virtual double finv(double w) const = 0;
    void initialize();  // Refactor common initialization logic into a protected method

    // Additional members and helper functions if needed
};

class LinearGrid : public Grid {
public:
    LinearGrid(double x0, double x1, int n);
private:
    double f(double x) const override {return  x; };
    double finv(double w) const override {return w; };
};

class QuadraticGrid : public Grid {
public:
    QuadraticGrid(double x0, double x1, int n);
private:
    double f(double x) const override;
    double finv(double w) const override;
};

class GeometricGrid : public Grid {
public:
    GeometricGrid(double x0, double x1, int n, double r = 1.1);
private:
    double r, logr, rn_minus_1;
    double f(double x) const override;
    double finv(double fx) const override;
};

class LogGrid : public Grid {
public:
    LogGrid(double x0, double x1, int n, double s = 10.0);
private:
    double s;
    double f(double x) const override;
    double finv(double fx) const override;
};

class AutomaticGrid : public Grid {
public:
    // AutomaticGrid(const Grid& grid, const std::vector<double>& f0, double delta = 1e-4);
    AutomaticGrid(const Grid& grid, const Eigen::VectorXd& f0, double delta = 1e-4);
private:
    double delta;
};

/*
// I've included a mkgrid function similar to the Python version, which will create grid objects based on the string kind.
// Grid* mkgrid(const std::string& kind, double x0, double x1, int n, ...); // Other parameters as needed
std::unique_ptr<Grid> mkgrid(const std::string& kind, ...); // Depending on how you plan to handle variable arguments
*/

}

#endif // GRID_H
