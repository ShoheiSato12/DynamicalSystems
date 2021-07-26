#include <fstream>



#ifndef RUNGEKUTTA4TH_H
#define RUNGEKUTTA4TH_H
#include<vector>
#include<cmath>

std::vector<long double> updateCoord(std::vector<long double>coord, std::vector<long double> increment, double step, double scale,int dimension);
std::vector<long double> rungeKutta4thSquare(std::vector<long double> (*function)(std::vector<long double>, double), 
                                  std::vector<long double>& coord, double param, double step, 
                                  int dimension );
void completeRungeKuttaToFile(std::vector<long double> (*function)(std::vector<long double>, double), std::vector<long double> initialCond,
                                                   double param, double step, int dimension, double time_span[2]);
std::vector<long double> rungeKutta4thSquarePertubation(std::vector<long double> (*function)(std::vector<long double>,std::vector<long double>, double), 
                                  std::vector<long double> coord,std::vector<long double> pertubation, double param, double step, 
                                  int dimension, std::vector<long double>& functionAval );
#endif