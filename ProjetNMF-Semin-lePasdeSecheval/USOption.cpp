//
//  USOption.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "USOption.hpp"

USOption::USOption(double strike, double maturity, payoffType payoff) : VanillaOption(strike, maturity, payoff) {}

bool USOption::isAmerican() const
{
    return true;
}

void testUSOption()
{
    double strike = 100;
    double maturity = 1;
    double stockPrice = 110;
    std::shared_ptr<VanillaOption> call1 = VanillaOption::newOption(strike, maturity, call, American);
    ASSERT_APPROX_EQUAL(call1->payoff(stockPrice), 10, 0.00001);
    
    stockPrice = 90;
    std::shared_ptr<VanillaOption> put1 = VanillaOption::newOption(strike, maturity, put, American);
    ASSERT_APPROX_EQUAL(put1->payoff(stockPrice), 10, 0.00001);
}
