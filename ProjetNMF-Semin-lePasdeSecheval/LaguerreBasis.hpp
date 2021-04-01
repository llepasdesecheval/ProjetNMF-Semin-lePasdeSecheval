//
//  LaguerreBasis.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined LAGUERREBASIS_HPP
#define LAGUERREBASIS_HPP

#include "stdafx.hpp"
#include "Basis.hpp"
#include "matlib.hpp"

class LaguerreBasis : public Basis
{
public:
    LaguerreBasis(std::size_t nbRegressors);
    
    double getKthRegressor(double x, std::size_t k) const override;
};

void testLaguerreBasis();

#endif /* LaguerreBasis_hpp */
