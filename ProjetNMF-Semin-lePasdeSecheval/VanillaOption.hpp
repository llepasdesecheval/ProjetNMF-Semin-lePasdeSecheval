//
//  VanillaOption.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined VANILLAOPTION_HPP
#define VANILLAOPTION_HPP

#include "stdafx.hpp"

// Useful enumerations
enum payoffType
{
    put = -1,
    call = 1
};

enum optionType
{
    American,
    European
};

class VanillaOption
{
public:
    /**
     Creates a vanilla option
     @param strike Strike of the option
     @param maturity Maturity of the option, in years
     @param payoff Call or put. Use the payoffType enumeration
     @param type American or European. Use the optionType enumeration.
     */
    static std::shared_ptr<VanillaOption> newOption(double strike, double maturity, payoffType payoff, optionType type);
    
    /**
     Virtual destructor
     */
    ~VanillaOption() {}
    
    double Strike() const;
    
    double Maturity() const;
    
    /**
     Returns the payoff of the option for a given stock price.
     @param stockPrice Stock price the payoff is applied to
     */
    double payoff(double stockPrice) const;
    
    /**
     Returns the vector of option payoff for a given stock price vector.
     @param stockPrices Stock price vector the payoff is applied to
     */
    std::vector<double> payoff(const std::vector<double>& stockPrices) const;
    
    virtual bool isAmerican() const = 0;
    
protected:
    VanillaOption(double strike, double maturity, payoffType payoff);
    double m_strike;
    double m_maturity;
    payoffType m_payoffType;
    optionType m_optionType;
};

void testVanillaOption();

#endif /* VanillaOption_hpp */
