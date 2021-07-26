
#include "../include/LyapExp.hpp"
#include "../include/LinearAlgebra.hpp"
#include <iostream>
std::vector<long double> lyapunovSpectrum(std::vector<long double> (*function)(std::vector<long double>, double),
                                std::vector<std::vector<long double>> (*jacobian)(std::vector<long double>&,double), 
                                std::vector<long double>& initialCond, double step, double param)
{
    // double step = 0.001;
    uint iterations = (uint)(100.0/step);
    std::vector<long double> coord = initialCond;
    std::vector<long double> lyapunovExponents (initialCond.size(),0);

    std::vector<std::vector<long double>> I=identityMatrix(initialCond.size());
    std::vector<std::vector<long double>> w=I;  
    std::vector<std::vector<long double>> J;
    std::vector<std::vector<long double>> auxJ(initialCond.size(),std::vector<long double> (initialCond.size(),0));
    for(uint i = 0; i < 100; i++)
    {
        coord = rungeKutta4thSquare(function, coord, param, step, initialCond.size());
    }

    for(uint i = 0; i < iterations; i++)
    {
        coord = rungeKutta4thSquare(function, coord, param, step, initialCond.size());
        auxJ=scalarXmat(step,jacobian(coord,param));
        J = matSum(I,auxJ);
        w = matMult(J,w);
        transpostSquare(w);
        gramSchmidt(w);
    
        for(uint j = 0; j < initialCond.size(); j++)
        {
            lyapunovExponents[j] +=  log(normOf(w[j]));
        }
    
        gramSchmidtNormal(w);
        transpostSquare(w);
    }
    for(uint j = 0; j < initialCond.size(); j++)
    {
        lyapunovExponents[j] = lyapunovExponents[j]/((long double)(100.0));
    }
    
    return lyapunovExponents;
}