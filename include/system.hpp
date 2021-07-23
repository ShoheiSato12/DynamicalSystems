#ifndef SYSTEMS_HPP
#define SYSTEMS_HPP
#include <vector>
#include "../include/LinearAlgebra.hpp"
std::vector<std::vector<double>> lorenzJacobian(std::vector<double>& coord, double rho);
std::vector<double> lorenz(std::vector<double> coord, double rho);
std::vector<double> classicalPendulum(std::vector<double> coord, double k);
std::vector<std::vector<double>> classicalPendulumJacobian(std::vector<double> &coord, double k);
std::vector<double> quantumPendulum(std::vector<double> coord, double gamma);
std::vector<std::vector<double>> quantumPendulumJacobian(std::vector<double> &coord, double gamma);
double fluctuation(double x, double y, double gamma);

#endif
