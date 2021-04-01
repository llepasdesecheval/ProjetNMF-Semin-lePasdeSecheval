//
//  USOption.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined USOPTION_HPP
#define USOPTION_HPP

#include "VanillaOption.hpp"

class USOption : public VanillaOption
{
public:
    /**
     Constructor for American option
     @param strike Option's strike price
     @param maturity Option's maturity
     @param payoff Option's payoff type; use payoffType enum
     */
    USOption(double strike, double maturity, payoffType payoff);
    
    bool isAmerican() const override;
};

void testUSOption();

#endif /* USOption_hpp */
