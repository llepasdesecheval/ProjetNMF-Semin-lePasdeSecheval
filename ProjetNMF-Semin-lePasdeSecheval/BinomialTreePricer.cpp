//
//  BinomialTreePricer.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "BinomialTreePricer.hpp"
#include "EurOption.hpp"

BinomialTreePricer::BinomialTreePricer(int nSteps) : m_nSteps(nSteps) {}

PricerOutput BinomialTreePricer::price(const VanillaOption& option, const BSMModel& model) const
{
    double step = (option.Maturity() - model.Date())/m_nSteps;
    double spot = model.StockPrice();
    double volatility = model.Volatility();
    double rate = model.InterestRate();
    double dividend = model.DividendYield();
    
    double upMove = exp(volatility * std::sqrt(step));
    double downMove = 1/upMove;
    double upProbability = (exp((rate - dividend) * step) - downMove)/(upMove - downMove);
    double stepDiscount = exp(-rate * step);
    double divCountBack = exp(dividend * step);
    
    // Terminal stock prices vector ; first element is highest terminal stock price
    std::vector<double> stockPrices(m_nSteps + 1, 0);
    // Option payoff vector
    std::vector<double> optionPayoff(m_nSteps + 1, 0);
    
    for (std::size_t i = 0; i <= m_nSteps; ++i)
    {
        double stockPrice = spot * std::pow(upMove, m_nSteps - i) * std::pow(downMove, i);
        stockPrices[i] = stockPrice;
        optionPayoff[i] = option.payoff(stockPrice);
    }
    
    for (std::size_t i = m_nSteps; i > 0; --i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            optionPayoff[j] = stepDiscount * (upProbability * (optionPayoff[j] - optionPayoff[j + 1]) + optionPayoff[j + 1]);
            
            if (option.isAmerican())
            {
                double stockPrice = stepDiscount * divCountBack * (upProbability * (stockPrices[j] - stockPrices[j + 1]) + stockPrices[j + 1]);
                stockPrices[j] = stockPrice;
                
                if (option.payoff(stockPrice) > optionPayoff[j])
                {
                    optionPayoff[j] = option.payoff(stockPrice);
                }
            }
        }
    }
    PricerOutput result(optionPayoff[0]);
    return result;
}

static void testBinomialTrees()
{
    BinomialTreePricer simpleTree(1);
    BSMModel model1(100, 0.1);
    std::shared_ptr<VanillaOption> EurCall = VanillaOption::newOption(100, 1, call, European);
    std::shared_ptr<VanillaOption> USCall = VanillaOption::newOption(100, 1, call, American);
    
    // On such a simple case, the two prices should be equal
    //ASSERT_APPROX_EQUAL(simpleTree.price(*EurCall, model1), simpleTree.price(*USCall, model1), 1E-6);
    ASSERT(simpleTree.price(*EurCall, model1).Estimate() == simpleTree.price(*USCall, model1).Estimate());
    ASSERT_APPROX_EQUAL(simpleTree.price(*USCall, model1).Estimate(), 4.995837, 1E-6);
    
    BinomialTreePricer largerTree(500);
    EurOption EurCall2(100, 1, call);
    /*
     Option on a non dividend-paying stock: US and Eur option prices should
     be roughly equal, and roughly equal to the Black-Scholes-Merton formula
     */
    ASSERT_APPROX_EQUAL(largerTree.price(*EurCall, model1).Estimate(), largerTree.price(*USCall, model1).Estimate(), 1E-6);
    ASSERT_APPROX_EQUAL(largerTree.price(*USCall, model1).Estimate(), EurCall2.price(model1).Estimate(), 0.01);
    
    // We now introduce dividends
    BSMModel model2(100, 0.1, 0, 0.03);
    ASSERT_APPROX_EQUAL(largerTree.price(*EurCall, model2).Estimate(), EurCall2.price(model2).Estimate(), 0.01);
    ASSERT_APPROX_EQUAL(largerTree.price(*USCall, model2).Estimate(), 2.925158, 1E-6);
}

void testBinomialTreePricer()
{
    TEST(testBinomialTrees);
}
