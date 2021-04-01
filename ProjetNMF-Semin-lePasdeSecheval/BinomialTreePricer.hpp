//
//  BinomialTreePricer.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined BINOMIALTREEPRICER_HPP
#define BINOMIALTREEPRICER_HPP

#include "Pricer.hpp"

/**
 Cox-Ross-Rubinstein binomial tree pricer
 */
class BinomialTreePricer : public Pricer
{
public:
    BinomialTreePricer(int nSteps);
    
    PricerOutput price(const VanillaOption& option, const BSMModel& model) const override;
    
private:
    int m_nSteps;
};

void testBinomialTreePricer();

#endif /* BinomialTreePricer_hpp */
