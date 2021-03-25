//
//  Basis.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 25/03/2021.
//

#if !defined BASIS_HPP
#define BASIS_HPP

#include "stdafx.hpp"
#include "Matrix.hpp"

/**
 Abstract class representing the regressor basis
 */
class Basis
{
public:
    /**
     Constructor
     @param nbRegressors Number of basis functions, excluding constant
     */
    Basis(std::size_t nbRegressors);
    
    /**
     Evaluates a basis function at the desired point
     @param x Evaluation point
     @param k Index of the basis function
     */
    virtual double getKthRegressor(double x, std::size_t k) const = 0;
    
    /**
     Returns the regressors matrix
     @param X Input vector
     */
    Matrix RegressorMatrix(const std::vector<double>& X) const;
    
    /**
     Returns the basis functions evaluated at the desired point
     @param x Evaluation point
     */
    std::vector<double> RegressorVector(double x) const;
    
protected:
    std::size_t m_nbRegressors;
};

#endif /* Basis_hpp */
