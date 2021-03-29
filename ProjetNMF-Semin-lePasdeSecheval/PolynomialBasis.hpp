//
//  PolynomialBasis.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 29/03/2021.
//

#if !defined POLYNOMIALBASIS_HPP
#define POLYNOMIALBASIS_HPP

#include "stdafx.hpp"
#include "Basis.hpp"

class PolynomialBasis : public Basis
{
public:
    PolynomialBasis(std::size_t nbRegressors);
    
    double getKthRegressor(double x, std::size_t k) const override;
};

void testPolynomialBasis();

#endif /* PolynomialBasis_hpp */
