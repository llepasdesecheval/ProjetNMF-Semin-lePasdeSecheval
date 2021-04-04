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
    Matrix result(m_nScenarii, m_nSteps + 1, false);

    // We compute all necessary parameters
    double drift = model.InterestRate() - model.DividendYield();
    double volatility = model.Volatility();
    double logDrift = (drift - 0.5 * volatility * volatility) * timeStep; // drift term of log price
    double logDiffusion = volatility * sqrt(timeStep); // diffusion term of log price
    result.setColumn(0, std::vector<double>(m_nScenarii, model.StockPrice()));

    for (std::size_t j = 0; j < m_nSteps; ++j)
    {
        std::vector<double> temp(result.ExtractColumn(j));
        // Log returns
        std::vector<double> returns = randN(m_nScenarii, logDiffusion, logDrift);
        // Returns
        std::transform(returns.begin(), returns.end(), returns.begin(), [](double x) { return exp(x); });
        // Price at step j + 1 = Price at step j * return(j, j + 1)
        std::transform(temp.begin(), temp.end(), returns.begin(), temp.begin(), std::multiplies<>{});
        result.setColumn(j + 1, temp);
    }
    return result;
}

Matrix MonteCarloPricer::LSWeights(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const
{
    double timeStep = (option.Maturity() - model.Date())/m_nSteps;
    double stepDiscount = exp(-model.InterestRate() * timeStep);
    double discount = 1; // Will change regularly
    
    Matrix result(m_basis->NbRegressors() + 1, m_nSteps - 1, false);
    
    // Computation of terminal payoffs
    Matrix payoffs(m_nScenarii, m_nSteps, true);
    payoffs.setColumn(m_nSteps - 1, option.payoff(pricePaths.ExtractColumn(m_nSteps)));
    
    // Backward induction and regressions
    for (std::size_t k = m_nSteps - 1; k > 0; --k)
    {
        std::vector<double> ITMStockPrices;
        std::vector<double> ITMFutureCashFlows;
        std::vector<std::size_t> ITMIndices;
        for (std::size_t i = 0; i < m_nScenarii; ++i)
        {
            double currentPrice = pricePaths(i, k);
            if (option.payoff(currentPrice) > 0) // Separate ITM paths
            {
                ITMIndices.push_back(i);
                ITMStockPrices.push_back(currentPrice);
                std::vector<double> temp = payoffs.ExtractRow(i);
                discount = 1;
                for (std::size_t m = 0; m < m_nSteps - k; ++m)
                {
                    discount *= stepDiscount;
                    temp[k + m] *= discount;
                }
                ITMFutureCashFlows.push_back(sum(temp));
            }
        }
        
        result.setColumn(k - 1, RegressionCoefficients(ITMStockPrices, ITMFutureCashFlows, *m_basis));
        std::vector<double> currentPrices = pricePaths.ExtractColumn(k);
        std::vector<double> continuation = m_basis->RegressorMatrix(currentPrices) * result.ExtractColumn(k - 1);
        ASSERT(continuation.size() == m_nScenarii);
        
        for (auto i : ITMIndices)
        {
            double exercise = option.payoff(currentPrices[i]);
            if (continuation[i] < exercise)
            {
                // There can be at most one exercise per path
                payoffs.setRow(i, std::vector<double>(m_nSteps, 0));
                payoffs(i, k - 1) = exercise;
            }
        }
    }
    return result;
}

PricerOutput MonteCarloPricer::LSPrice(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const
{
    ASSERT(m_basis != nullptr);
    
    double timeStep = (option.Maturity() - model.Date())/m_nSteps;
    double stepDiscount = exp(-model.InterestRate() * timeStep);
    double discount = 1; // Will be changed regularly

    // Step 1 : Computation of terminal payoffs
    Matrix payoffs(m_nScenarii, m_nSteps, true);
    payoffs.setColumn(m_nSteps - 1, option.payoff(pricePaths.ExtractColumn(m_nSteps)));

    // Step 2 : Backward induction and regressions
    for (std::size_t k = m_nSteps - 1; k > 0; --k)
    {
        std::vector<double> ITMStockPrices;
        std::vector<double> ITMFutureCashFlows;
        std::vector<std::size_t> ITMIndices;
        for (std::size_t i = 0; i < m_nScenarii; ++i)
        {
            double currentPrice = pricePaths(i, k);
            if (option.payoff(currentPrice) > 0) // Separate ITM paths
            {
                ITMIndices.push_back(i);
                ITMStockPrices.push_back(currentPrice);
                std::vector<double> temp = payoffs.ExtractRow(i);
                discount = 1;
                for (std::size_t m = 0; m < m_nSteps - k; ++m)
                {
                    discount *= stepDiscount;
                    temp[k + m] *= discount;
                }
                ITMFutureCashFlows.push_back(sum(temp));
            }
        }
        
        std::vector<double> regCoef = RegressionCoefficients(ITMStockPrices, ITMFutureCashFlows, *m_basis);
        std::vector<double> currentPrices = pricePaths.ExtractColumn(k);
        std::vector<double> continuation = m_basis->RegressorMatrix(currentPrices) * regCoef;
        ASSERT(continuation.size() == m_nScenarii);
        
        for (auto i : ITMIndices)
        {
            double exercise = option.payoff(currentPrices[i]);
            if (continuation[i] < exercise)
            {
                // There can be at most one exercise per path
                payoffs.setRow(i, std::vector<double>(m_nSteps, 0));
                payoffs(i, k - 1) = exercise;
            }
        }
    }

    // Step 3 : Discounting option's cash flows and averaging across scenarios
    std::vector<double> USOptionDCF;
    for (std::size_t i = 0; i < m_nScenarii; ++i)
    {
        discount = 1;
        for (std::size_t j = 0; j < m_nSteps; ++j)
        {
            // We discount the cash flow by the appropriate amount
            discount *= stepDiscount;
            payoffs(i, j) *= discount;
        }
        // We add it to our payoff sample vector
        USOptionDCF.push_back(max(payoffs.ExtractRow(i)));
    }
    
    // We compare it to exercising at time 0
    std::transform(USOptionDCF.begin(), USOptionDCF.end(), USOptionDCF.begin(), [&model, &option](double x)
                   { return (option.payoff(model.StockPrice()) > x) ? option.payoff(model.StockPrice()) : x; });

    PricerOutput result(USOptionDCF, 0.99);
    
    return result;
}

PricerOutput MonteCarloPricer::ABPrice(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const
{
    return PricerOutput(0);
}

PricerOutput MonteCarloPricer::price(const VanillaOption& option, const BSMModel& model) const
{
    if (option.isAmerican())
    {
        // Step 1 : Generation of stock price paths
        Matrix pricePaths = generatePricePaths(model, option.Maturity());
        
        return LSPrice(option, model, pricePaths);
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

static void testUSOptionPricing()
{
     LaguerreBasis basis(2);
     MonteCarloPricer mc(25000, 250, &basis);
     BSMModel model1(100, 0.2, 0.05, 0.08); // High dividend yield so the call might be exercised
     std::shared_ptr<VanillaOption> USCall = VanillaOption::newOption(100, 1, call, American);
     BinomialTreePricer tree(52);
     
     PricerOutput outputCall = mc.price(*USCall, model1);
     double comparisonCall = tree.price(*USCall, model1).Estimate();
     ASSERT_APPROX_EQUAL(outputCall.Estimate(), comparisonCall, outputCall.Confidence());
     
     BSMModel model2(100, 0.2, 0.05); // No dividends, so the put is more likely to be exercised
     std::shared_ptr<VanillaOption> USPut = VanillaOption::newOption(100, 1, put, American);
     PricerOutput outputPut = mc.price(*USPut, model2);
     double comparisonPut = tree.price(*USPut, model2).Estimate();
     ASSERT_APPROX_EQUAL(outputPut.Estimate(), comparisonPut, outputPut.Confidence());
}

void testMonteCarloPricer()
{
    TEST(testGeneratePath);
    TEST(testEuropeanOptionPricing);
    TEST(testUSOptionPricing);
}
