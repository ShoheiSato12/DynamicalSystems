#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP
#include <vector>
#include <cmath>

//Very basic operations
inline long double dotProduct(std::vector<long double> &v1, std::vector<long double> &v2);
long double normOf(std::vector<long double> &vec);
inline double normSquare(std::vector<long double> &vec);
inline void normalize(std::vector<long double> &vec);

//Operations on matrices
void gramSchmidt(std::vector<std::vector<long double>>& matrix);
void gramSchmidtNormal(std::vector<std::vector<long double>>& matrix);
void printMatrix(std::vector<std::vector<long double>>& matrix);
void transpostSquare(std::vector<std::vector<long double>>& matrix);

//auxiliary operations on matrices
std::vector<long double> projectionIntoU(std::vector<long double> vectorV,std::vector<long double> vectorU);
std::vector<long double>  sumOfProjections(uint position, std::vector<std::vector<long double>>& matrix);

std::vector<std::vector<long double>> matMult(std::vector<std::vector<long double>>& leftMatrix, std::vector<std::vector<long double>>& rigthMatrix);
std::vector<std::vector<long double>> identityMatrix(uint order);
std::vector<std::vector<long double>> matSum(std::vector<std::vector<long double>>& leftMatrix, std::vector<std::vector<long double>>& rigthMatrix);
std::vector<std::vector<long double>> scalarXmat(double num,std::vector<std::vector<long double>> matrix);

#endif