#include "../include/rungekutta4thSquare.hpp"


std::vector<long double> rungeKutta4thSquare(std::vector<long double> (*function)(std::vector<long double>, double), 
                                  std::vector<long double> &coord, double param, double step, 
                                  int dimension )
{
    std::vector<long double> k1 (dimension);
    std::vector<long double> k2 (dimension);
    std::vector<long double> k3 (dimension);
    std::vector<long double> k4 (dimension);
    
    k1 = function(coord, param);
    k2 = function(updateCoord(coord, k1, (double)step, 2.0,dimension), param);
    k3 = function(updateCoord(coord, k2, (double)step, 2.0,dimension), param);
    k4 = function(updateCoord(coord, k3, (double)step, 1.0,dimension), param);
    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + ((double)step/(double)6.0)*(k1[i]+2.0*k2[i] + 2.0*k3[i]+ k4[i]);
    }
    return coord;
}
std::vector<long double> rungeKutta4thSquarePertubation(std::vector<long double> (*function)(std::vector<long double>,std::vector<long double>, double), 
                                  std::vector<long double> coord,std::vector<long double> pertubation, double param, double step, 
                                  int dimension, std::vector<long double>& functionAval )
{
    std::vector<long double> k1 (dimension);
    std::vector<long double> k2 (dimension);
    std::vector<long double> k3 (dimension);
    std::vector<long double> k4 (dimension);
    
    k1 = function(coord, pertubation, param);
    k2 = function(updateCoord(coord, k1, (double)step, 2.0,dimension), pertubation, param);
    k3 = function(updateCoord(coord, k2, (double)step, 2.0,dimension), pertubation, param);
    k4 = function(updateCoord(coord, k3, (double)step, 1.0,dimension), pertubation, param);
    for(int i = 0; i < dimension; i++)
    {
        coord[i] = coord[i] + ((double)step/(double)6.0)*(k1[i]+2.0*k2[i] + 2.0*k3[i]+ k4[i]);
    }
    functionAval = k4;
    return coord;
}
std::vector<long double> updateCoord(std::vector<long double> coord, std::vector<long double> increment, double step, double scale,int dimension)
{
    
    std::vector<long double> updatedCoord (dimension);
    for(int i = 0; i < dimension; i++)
    {
        updatedCoord[i] = coord[i] + ((double)step/(double)scale)*(increment[i]);
    }
    return updatedCoord;
}

void completeRungeKuttaToFile(std::vector<long double> (*function)(std::vector<long double>, double), std::vector<long double> initialCond,
                                                   double param, double step, int dimension, double time_span[2])
{
    int iterations = (int)(fabs(time_span[1]-time_span[0])/step);
    std::ofstream fileName;
    fileName.open("ouputrk4th.dat");
    //std::vector<long double> auxVec = rungeKutta4thSquare(function, initialCond, param, step, dimension);
    for(int j = 0; j < iterations; j++)
    {
        initialCond = rungeKutta4thSquare(function, initialCond, param, step, dimension);
        for(int i = 0; i < dimension; i++)
        {
            fileName << initialCond[i] <<"  ";
        }
        fileName<< std::endl;
    }
    fileName.close();

}