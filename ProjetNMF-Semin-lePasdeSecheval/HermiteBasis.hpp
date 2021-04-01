//
//  HermiteBasis.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined HERMITEBASIS_HPP
#define HERMITEBASIS_HPP

#include "stdafx.hpp"
#include "Basis.hpp"
#include "matlib.hpp"

class HermiteBasis : public Basis
{
public:
    HermiteBasis(std::size_t nbRegressors);
    
    double getKthRegressor(double x, std::size_t k) const override;
};

void testHermiteBasis();

#endif /* HermiteBasis_hpp */
