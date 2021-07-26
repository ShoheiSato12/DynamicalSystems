#include "../include/LinearAlgebra.hpp"
#include<iostream>
inline long double dotProduct(std::vector<long double> &v1, std::vector<long double> &v2)
{
    long double result =  0;
    for(uint i = 0; i < v1.size(); i++)
    {
        result += v1[i]*v2[i];
    }
    return result;
}
long double normOf(std::vector<long double> &vec)
{
    return sqrt(dotProduct(vec,vec));
}
inline double normSquare(std::vector<long double> &vec)
{
    return dotProduct(vec,vec);
}
inline void normalize(std::vector<long double> &vec)
{
    double norm = normOf(vec);
    for(uint i = 0; i < vec.size(); i ++)
    {
        vec[i]/=norm;
    }
}
std::vector<long double> projectionIntoU(std::vector<long double> vectorV,std::vector<long double> vectorU)
{
    double coef =  dotProduct(vectorU,vectorV)*(1.0/normSquare(vectorU));
    std::vector<long double> proj (vectorU.size());
    for(uint i = 0; i < vectorU.size(); i++)
    {
        proj[i] = coef * vectorU[i];
    }
    return proj;
}
 std::vector<long double>  sumOfProjections(uint position, std::vector<std::vector<long double>>& matrix)
{
    std::vector<long double>  sum(matrix[0].size(),0);
    
    for(uint i = 0; i < position; i++)
    {
   

        for(uint j = 0; j < matrix[0].size(); j++)
        {
            sum[j] -=  projectionIntoU(matrix[position], matrix[i])[j];
        }

    }
    for(uint j = 0; j < matrix[0].size(); j++)
    {
        sum[j] = matrix[position][j] + sum[j];
    }
    return sum;
}
void transpostSquare(std::vector<std::vector<long double>>& matrix)
{
    std::vector<std::vector<long double>> aux = matrix; 
    for(uint i = 0; i < matrix.size(); i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
             matrix[i][j] = aux[j][i];
        }
    }
  
}
void gramSchmidt(std::vector<std::vector<long double>>& matrix)
{
    
    for(uint i = 1; i < matrix[0].size(); i ++)
    {
        matrix[i] = sumOfProjections(i, matrix);     
    }
}
void gramSchmidtNormal(std::vector<std::vector<long double>>& matrix)
{
    gramSchmidt(matrix);
    std::vector<long double> norms (matrix.size());
    for(uint j = 0; j < matrix.size(); j++)
    {
        norms[j] = normOf(matrix[j]);
    }
    for(uint i = 0; i < matrix[0].size(); i++)
    {
        for(uint j = 0; j < matrix.size(); j++)
        {
            matrix[j][i] = matrix[j][i]/norms[j];

        }
    }
}
void printMatrix(std::vector<std::vector<long double>>& matrix)
{
    for(uint i = 0; i < matrix.size(); i++)
    {
        for(uint j = 0; j < matrix[0].size(); j++)
        {
            std::cout<<matrix[i][j]<<" ";
        }
        std::cout<<"\n";
    }
}
std::vector<std::vector<long double>> identityMatrix(uint order)
{
    std::vector<std::vector<long double>> id (order, std::vector<long double>(order,0));
    for(uint i = 0; i < order; i++)
    {
        id[i][i] = 1.0;
    }
    return id;
}  
std::vector<std::vector<long double>> matMult(std::vector<std::vector<long double>>& leftMatrix, std::vector<std::vector<long double>>& rigthMatrix) 
{
    std::vector<std::vector<long double>> result (leftMatrix.size(),std::vector<long double>(rigthMatrix[0].size(),0));

   
    for (uint i = 0; i < leftMatrix.size(); i++) {
        for (uint j = 0; j < rigthMatrix[0].size(); j++) {
            result[i][j] = 0;
            for (uint k = 0; k < rigthMatrix.size(); k++)
                result[i][j] += leftMatrix[i][k] * rigthMatrix[k][j];
        }
    }
    return result;
 
}
std::vector<std::vector<long double>> matSum(std::vector<std::vector<long double>>& leftMatrix, std::vector<std::vector<long double>>& rigthMatrix)
{
    std::vector<std::vector<long double>> result (leftMatrix.size(),std::vector<long double> (leftMatrix.size(),0));
    for(uint i=0; i<leftMatrix.size();i++)
    {
        for(uint j=0;j<leftMatrix[0].size();j++)
        {
            result[i][j]=leftMatrix[i][j]+rigthMatrix[i][j];
        }
    }
    return result;
}
std::vector<std::vector<long double>> scalarXmat(double num,std::vector<std::vector<long double>> matrix)
{
    std::vector<std::vector<long double>> result(matrix.size(),std::vector<long double>(matrix[0].size(),0));
    for(uint i=0;i<matrix.size();i++)
    {
        for(uint j=0;j<matrix[0].size();j++)
        {
            result[i][j]=num*matrix[i][j];
        }
    }
    return result;
}

