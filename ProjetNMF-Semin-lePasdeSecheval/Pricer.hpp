//
//  Pricer.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 29/03/2021.
//

#if !defined PRICER_HPP
#define PRICER_HPP

#include "stdafx.hpp"
#include "PricerOutput.hpp"
#include "VanillaOption.hpp"
#include "BSMModel.hpp"
#include "PricerOutput.hpp"

/**
 Interface class for vanilla option pricers
 */
class Pricer
{
public:
    /** Virtual destructor */
    ~Pricer() {}
    
    /**
     Pricing method
     @param option Vanilla option to price
     @param model Black-Scholes model containing the pricing parameters
     */
    virtual PricerOutput price(const VanillaOption& option, const BSMModel& model) const = 0;
};

#endif /* Pricer_h */
