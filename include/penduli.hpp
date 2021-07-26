#include <vector>
#include <cstdlib>
#include<ctime>
#include <cmath>

#ifndef PENDULI_H
#define PENDULI_H

std::vector<long double> classicalPendulum(std::vector<long double> coord, double k);
std::vector<long double> quantumPendulum(std::vector<long double> coord, double gamma);
double flutuation(double x, double y, double gamma);
#endif