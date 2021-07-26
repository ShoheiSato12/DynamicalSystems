#ifndef LYAPUNOV_H
#define LYAPUNOV_V
#include "../include/LinearAlgebra.hpp"
#include "../include/rungekutta4thSquare.hpp"
#include <ctime>
#include <cstdlib>
#include<iostream>

std::vector<long double> lyapunovSpectrum(std::vector<long double> (*function)(std::vector<long double>, double),
                                std::vector<std::vector<long double>> (*jacobian)(std::vector<long double>&,double), 
                                std::vector<long double>& initalCond, double step, double param);
#endif
