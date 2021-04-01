//
//  EurOption.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "EurOption.hpp"

EurOption::EurOption(double strike, double maturity, payoffType payoff) : VanillaOption(strike, maturity, payoff) {}

bool EurOption::isAmerican() const
{
    return false;
}

PricerOutput EurOption::price(const BSMModel& model) const
{
    // Ex-dividend spot
    double timeToMaturity = m_maturity - model.Date();
    double forwardSpot = model.StockPrice() * exp(-model.DividendYield() * timeToMaturity);
    // Discount factor
    double discount = exp(-model.InterestRate() * timeToMaturity);
    
    double totalVol = model.Volatility() * sqrt(timeToMaturity);
    double d1 = log(forwardSpot/(m_strike * discount))/totalVol + 0.5 * totalVol;

    PricerOutput result(m_payoffType * (forwardSpot * NormSDist(m_payoffType * d1) - m_strike * discount * NormSDist(m_payoffType * (d1 - totalVol))));
    
    return result;
}

void testEurOption()
{
    EurOption callOption(105, 2, call);
    EurOption putOption(95, 2, put);
    
    BSMModel model(100, 0.25, 0.05, 0.02, 1.0);
    ASSERT_APPROX_EQUAL(callOption.price(model).Estimate(), 8.941176, 0.000001);
    ASSERT_APPROX_EQUAL(putOption.price(model).Estimate(), 6.031656, 0.000001);
}
