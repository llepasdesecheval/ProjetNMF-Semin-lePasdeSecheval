//
//  Basis.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 25/03/2021.
//

#include "Basis.hpp"

Basis::Basis(std::size_t nbRegressors) : m_nbRegressors(nbRegressors) {}

Matrix Basis::RegressorMatrix(const std::vector<double>& X) const
{
    std::size_t n = X.size();
    Matrix L(n, m_nbRegressors + 1, false);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < m_nbRegressors + 1; ++j)
        {
            L(i, j) = getKthRegressor(X[i], j);
        }
    }
    return L;
}

std::vector<double> Basis::RegressorVector(double x) const
{
    std::vector<double> result;
    for (std::size_t j = 0; j < m_nbRegressors + 1; ++j)
    {
        result.push_back(getKthRegressor(x, j));
    }
    return result;
}
