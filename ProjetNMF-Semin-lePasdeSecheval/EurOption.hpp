//
//  EurOption.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined EUROPTION_HPP
#define EUROPTION_HPP

#include "VanillaOption.hpp"
#include "BSMModel.hpp"
#include "PricerOutput.hpp"

class EurOption : public VanillaOption
{
public:
    /**
     Constructor for European option
     @param strike Option's strike price
     @param maturity Option's maturity
     @param payoff Option's payoff type; use payoffType enum
     */
    EurOption(double strike, double maturity, payoffType payoff);
    
    /**
     Test method: returns the option's analytic price in the Black-Scholes-Merton framework
     @param model Black-Scholes-Merton model used to price the option
     */
    PricerOutput price(const BSMModel& model) const;
    
    bool isAmerican() const override;
};

void testEurOption();

#endif /* EurOption_hpp */
