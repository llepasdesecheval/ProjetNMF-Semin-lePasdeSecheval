//
//  MonteCarloPricer.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "MonteCarloPricer.hpp"
#include "BinomialTreePricer.hpp"
#include "EurOption.hpp"
#include "USOption.hpp"
#include "PolynomialBasis.hpp"
#include "HermiteBasis.hpp"
#include "LaguerreBasis.hpp"

MonteCarloPricer::MonteCarloPricer(std::size_t nScenarii, std::size_t nSteps, Basis* basis) :
m_nScenarii(nScenarii), m_nSteps(nSteps), m_basis(basis) {}

Matrix MonteCarloPricer::generatePricePaths(const BSMModel& model, double toDate) const
{
    double timeStep = (toDate - model.Date())/m_nSteps;
    Matrix res(m_nScenarii, m_nSteps + 1, false);

    // We compute all necessary parameters
    double drift = model.InterestRate() - model.DividendYield();
    double volatility = model.Volatility();
    double logDrift = (drift - 0.5 * volatility * volatility) * timeStep; // drift term of log price
    double logDiffusion = volatility * sqrt(timeStep); // diffusion term of log price
    res.setColumn(0, std::vector<double>(m_nScenarii, model.StockPrice()));

    for (std::size_t j = 0; j < m_nSteps; ++j)
    {
        std::vector<double> temp(res.ExtractColumn(j));
        // Log returns
        std::vector<double> returns = randN(m_nScenarii, logDiffusion, logDrift);
        // Returns
        std::transform(returns.begin(), returns.end(), returns.begin(), [](double x) { return exp(x); });
        // Price at step j + 1 = Price at step j * return(j, j + 1)
        std::transform(temp.begin(), temp.end(), returns.begin(), temp.begin(), std::multiplies<>{});
        res.setColumn(j + 1, temp);
    }
    return res;
}

PricerOutput MonteCarloPricer::price(const VanillaOption& option, const BSMModel& model) const
{
    if (option.isAmerican())
    {
        PricerOutput result(0);

        return result;
    }
    else
    {
        // To price a European option, we only need the spot at maturity
        MonteCarloPricer replace(m_nScenarii, 1);
        Matrix priceSample = replace.generatePricePaths(model, option.Maturity());
        std::vector<double> payoffs = option.payoff(priceSample.ExtractColumn(1));
        double discount = exp(-model.InterestRate() * (option.Maturity() - model.Date()));
        std::transform(payoffs.begin(), payoffs.end(), payoffs.begin(), [discount](double x) { return discount * x; });
        PricerOutput result(payoffs, 0.99);
        return result;
    }
}

static void testGeneratePath() {

    BSMModel model(100, 0.1, 0.05, 0, 1, 0); // Low volatility to have narrow final prices sample

    std::size_t nPaths = 10000;
    std::size_t nSteps = 5;
    MonteCarloPricer mc(nPaths, nSteps);
    double maturity = 2.0;

    Matrix pricePaths = mc.generatePricePaths(model, maturity);
    std::vector<double> finalPrices = pricePaths.ExtractColumn(5);
    ASSERT_APPROX_EQUAL(mean(finalPrices), exp(model.InterestRate() * (maturity - model.Date())) * model.StockPrice(), 0.5);
}

static void testEuropeanOptionPricing()
{
    MonteCarloPricer mc(1000000, 10);
    BSMModel model1(100, 0.1);
    std::shared_ptr<VanillaOption> EurCall = VanillaOption::newOption(100, 1, call, European);
    EurOption EurCall2 = EurOption(100, 1, call);
    PricerOutput output = mc.price(*EurCall, model1);
    ASSERT_APPROX_EQUAL(output.Estimate(), EurCall2.price(model1).Estimate(), output.Confidence());
}

void testMonteCarloPricer()
{
    TEST(testGeneratePath);
    TEST(testEuropeanOptionPricing);
}
