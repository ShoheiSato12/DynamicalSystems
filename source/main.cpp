#include<iostream>
#include<fstream>
#include<iterator> // for iterators
#include<vector> // for vectors
#include<cmath>
// #include "../include/gnuplot-iostream.h"
 //#include "../include/rungekutta4thSquare.hpp"
//#include "../include/biffurcation.hpp"
//#include "../include/penduli.hpp"
#include "../include/LyapExp.hpp"
#include "../include/system.hpp"

int main()
{
    std::vector<double> integrationAux = {6.0,1.0,0.0,1.0};
    double parameter = 0;
    for (int i = 0; i < 10; i++)
    {
        std::vector<long double> lya = lyapunovSpectrum(quantumPendulum, quantumPendulumJacobian, integrationAux, 0.001, parameter);
        // std::vector<std::vector<double>> A = matMult(M,N);
        std::cout << "--- Gamma = " << parameter << "---------" << std::endl;
        std::cout << "lyapunov exponents: " << (lya[0]) << ", " << (lya[1]) << ", " << (lya[2]) << "\n ";
        std::cout<<"lyapunov numbers: "<<exp(lya[0])<<", "<<exp(lya[1])<<", "<<exp(lya[2])<<"\n ";

        std::cout<<"sum of lyapunov exponents: "<< lya[0] + lya[1] + lya[2];
        parameter += 0.0001;
    }
    return 0;    
}
