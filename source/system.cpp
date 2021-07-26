#include<time.h>
#include "../include/system.hpp"
//Lorenz
std::vector<long double> lorenz(std::vector<long double> coord, double rho)
{
    std::vector<long double> xcoord (3,0);
    double a=10,b=8/3;
    xcoord[0] = a*coord[1] - 10.0 * coord[0];
    xcoord[1] = rho*coord[0] - coord[1] - coord[0]*coord[2];
    xcoord[2] = coord[0]*coord[1] -b*coord[2];
    return xcoord;
}
std::vector<std::vector<long double>> lorenzJacobian(std::vector<long double>& coord,double rho)
{
    std::vector<std::vector<long double>> jacobian (3,std::vector<long double>(3,0));
    double a=10,b=8/3;
    jacobian[0][0] = -a; 
    jacobian[0][1] =  a;
    jacobian[0][2] =  0.0;  
    jacobian[1][0] = -coord[2] + rho;
    jacobian[1][1] = -1;
    jacobian[1][2] = -coord[0];
    jacobian[2][0] = coord[1];
    jacobian[2][1] = coord[0];
    jacobian[2][2] = -b;
    return jacobian;
}
//Classical pendulum
std::vector<long double> classicalPendulum(std::vector<long double> coord, double k)
{
    std::vector<long double> coord_dot (4);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = -coord[0] - k*coord[0]*coord[1]*coord[1];
    coord_dot[3] = -coord[1] - k*coord[1]*coord[0]*coord[0];
    
    return coord_dot;
}
std::vector<std::vector<long double>> classicalPendulumJacobian(std::vector<long double> &coord, double k)
{
    std::vector<std::vector<long double>> jacobian (coord.size(),std::vector<long double>(coord.size(),0));
    jacobian[0] = {0, 0 , 1 , 0};
    jacobian[1] = {0 , 0 , 0 , 1};
    jacobian[2] = {-1-k*coord[1]*coord[1] , -2*k*coord[0]*coord[1] , 0 , 0};
    jacobian[3] = {-2*k*coord[0]*coord[1] , -1-k*coord[0]*coord[0] , 0 , 0};
    return jacobian;
}
//Quantum pendulum
std::vector<long double> quantumPendulum(std::vector<long double> coord, double gamma)
{
    std::vector<long double> coord_dot (4,0);
    srand (time(NULL));
    double f = fluctuation(coord[0],coord[1], gamma);
    coord_dot[0] = coord[2];
    coord_dot[1] = coord[3];
    coord_dot[2] = -gamma*coord[2] - coord[0] - coord[0]*coord[1]*coord[1] + f*rand();
    coord_dot[3] = -gamma*coord[3] - coord[1] - coord[1]*coord[0]*coord[0] + f*rand();

    return coord_dot; 
}
std::vector<std::vector<long double>> quantumPendulumJacobian(std::vector<long double> &coord, double gamma)
{
    srand (time(NULL));
    double f = fluctuation(coord[0],coord[1], gamma);
    std::vector<std::vector<long double>> jacobian (coord.size(),std::vector<long double>(coord.size(),0));
    jacobian[0] = {0 , 0 , 1 , 0};
    jacobian[1] = {0 , 0 , 0 , 1};
    jacobian[2] = {-1-coord[1]*coord[1]-gamma*coord[2]+f*rand() , -2*coord[0]*coord[1]-gamma*coord[2]+f*rand() , 0 , 0};
    jacobian[3] = {-2*coord[0]*coord[1]-gamma*coord[3]+f*rand() , -1-coord[0]*coord[0]-gamma*coord[3]+f*rand() , 0 , 0};
    return jacobian;

}
double fluctuation(double x, double y, double gamma)
{
        return sqrt(
            gamma*(0.62832912000*(pow(x,2.0) + pow(y,2.0)) +
             0.01205153100* pow(x,2.0) * pow(y,2.0)  +
             0.56437351570*(pow(x,4.0) + pow(y,4.0)) +
             0.13998728990*(pow(x,6.0) + pow(y,6.0)) +
             0.06202680930*(pow(x,2.0)* pow(y,4.0) +  pow(y,2.0)*pow(x,4.0)) +
             0.03555913955* pow(x,4.0) *  pow(y,4.0)  +
             0.01777956970*( pow(x,8.0) +  pow(y,8.0)) + 1.155436517)/
             pow((1.0 + 0.3345167463*(pow(x,2.0) +  pow(y,2.0)) +
              0.09871707060*(pow(x,4.0) +  pow(y,4.0))),2.0));
}

