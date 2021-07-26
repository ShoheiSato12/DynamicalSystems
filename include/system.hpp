#ifndef SYSTEMS_HPP
#define SYSTEMS_HPP
#include <vector>
#include "../include/LinearAlgebra.hpp"
std::vector<std::vector<long double>> lorenzJacobian(std::vector<long double>& coord, double rho);
std::vector<long double> lorenz(std::vector<long double> coord, double rho);
std::vector<long double> classicalPendulum(std::vector<long double> coord, double k);
std::vector<std::vector<long double>> classicalPendulumJacobian(std::vector<long double> &coord, double k);
std::vector<long double> quantumPendulum(std::vector<long double> coord, double gamma);
std::vector<std::vector<long double>> quantumPendulumJacobian(std::vector<long double> &coord, double gamma);
double fluctuation(double x, double y, double gamma);

#endif
