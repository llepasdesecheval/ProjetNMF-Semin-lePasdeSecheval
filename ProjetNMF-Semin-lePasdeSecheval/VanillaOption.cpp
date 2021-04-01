//
//  VanillaOption.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "VanillaOption.hpp"
#include "EurOption.hpp"
#include "USOption.hpp"

VanillaOption::VanillaOption(double strike, double maturity, payoffType payoff) :
m_strike(strike), m_maturity(maturity), m_payoffType(payoff) {}

double VanillaOption::payoff(double stockPrice) const
{
    return (m_payoffType * stockPrice > m_payoffType * m_strike) ?
    m_payoffType * (stockPrice - m_strike) : 0;
}

std::vector<double> VanillaOption::payoff(const std::vector<double>& stockPrices) const
{
    std::vector<double> res(stockPrices.size());
    /*
    std::for_each(res.begin(), res.end(), [this](double x){ payoff(x); });
    double strike = m_strike;
    payoffType payoff = m_payoffType;
    std::for_each(res.begin(), res.end(), [strike, payoff](double x)
    { return (payoff * x > payoff * strike) ? payoff * (x - strike) : 0; });
     */
    for (std::size_t i = 0; i < stockPrices.size(); ++i)
    {
        res[i] = payoff(stockPrices[i]);
    }
    return res;
}

double VanillaOption::Maturity() const
{
    return m_maturity;
}

double VanillaOption::Strike() const
{
    return m_strike;
}

std::shared_ptr<VanillaOption> VanillaOption::newOption(double strike, double maturity, payoffType payoff, optionType type)
{
    std::shared_ptr<VanillaOption> ret;
    if (type == American)
    {
        ret = std::make_shared<USOption>(strike, maturity, payoff);
    }
    else
    {
        ret = std::make_shared<EurOption>(strike, maturity, payoff);
    }
    return ret;
}

static void testPayoff()
{
    std::shared_ptr<VanillaOption> EurCall = VanillaOption::newOption(100, 1, call, European);
    std::vector<double> v = {110, 102, 98, 103.5, 97.89, 111.234897};
    std::vector<double> expected = {10, 2, 0, 3.5, 0, 11.234897};
    std::vector<double> output = EurCall->payoff(v);
    for (std::size_t i = 0; i < v.size(); ++i)
    {
        ASSERT_APPROX_EQUAL(output[i], expected[i], 1e-6);
    }
}

void testVanillaOption()
{
    TEST(testEurOption);
    TEST(testUSOption);
    TEST(testPayoff);
}
