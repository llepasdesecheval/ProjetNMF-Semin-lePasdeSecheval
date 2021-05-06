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
#include "LineChart.hpp"
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

std::vector<double> MonteCarloPricer::nestedMC(double St, double horizon, std::size_t nSuccessors, const BSMModel& model) const
{
    double volatility = model.Volatility();
    std::vector<double> result = randN(nSuccessors, volatility * sqrt(horizon),
                                        (model.InterestRate() - model.DividendYield()
                                         - 0.5 * volatility * volatility) * horizon);
    // Returns
    std::transform(result.begin(), result.end(), result.begin(), [St](double x) { return St * exp(x); });
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
    double timeStep = (option.Maturity() - model.Date())/m_nSteps;
    double stepDiscount = exp(-model.InterestRate() * timeStep);
    
    Matrix M(pricePaths.nRows(), pricePaths.nCols(), true); // Martingale matrix
    Matrix U(pricePaths.nRows(), pricePaths.nCols(), true); // Upper bound matrix
    Matrix weights = LSWeights(option, model, pricePaths);
    for (std::size_t j = 1; j < m_nSteps; ++j)
    {
        for (std::size_t i = 0; i < m_nScenarii; ++i)
        {
            double spot = pricePaths(i, j);
            std::vector<double> stepWeights = weights.ExtractColumn(j - 1);
            double payoff = option.payoff(spot);
            double continuation = std::inner_product(stepWeights.begin(), stepWeights.end(), m_basis->RegressorVector(spot).begin(), 0);
            double Vt = (payoff < continuation) ? continuation : payoff;
            std::vector<double> nestedSim = nestedMC(spot, timeStep, 100, model);
            std::vector<double> nestedContinuations = m_basis->RegressorMatrix(nestedSim) * stepWeights;
            std::vector<double> nestedPayoffs(nestedSim);
            std::transform(nestedPayoffs.begin(), nestedPayoffs.end(), nestedPayoffs.begin(), [&option](double x)
                           { return option.payoff(x); });
            std::vector<double> temp(nestedPayoffs);
            std::transform(nestedPayoffs.begin(), nestedPayoffs.end(), nestedContinuations.begin(), temp.begin(), std::greater<double> {});
            double VtEst = mean(temp);
            M(i, j) = M(i, j - 1) / stepDiscount + (Vt - VtEst); // Optimal martingale
            if (j == m_nSteps)
            {
                U(i, j) = (U(i, j - 1)/ stepDiscount > mean(nestedPayoffs) - M(i, j))  ? U(i, j - 1)/ stepDiscount : mean(nestedPayoffs) - M(i, j);
            }
            else
            {
                U(i, j) = (U(i, j - 1)/ stepDiscount > option.payoff(spot) - M(i, j))  ? U(i, j - 1)/ stepDiscount : option.payoff(spot) - M(i, j);
            }
        }
    }
    
    for (std::size_t i = 0; i < m_nScenarii; ++i)
    {
        double spot = pricePaths(i, m_nSteps);
        double terminalPayoff = option.payoff(spot);
        M(i, m_nSteps) = M(i, m_nSteps - 1) / stepDiscount; // Optimal martingale
        U(i, m_nSteps) = (U(i, m_nSteps - 1)/ stepDiscount > terminalPayoff - M(i, m_nSteps)) ? U(i, m_nSteps - 1)/ stepDiscount : terminalPayoff - M(i, m_nSteps);
    }
    
    double discount = exp(-model.InterestRate() * option.Maturity());
    std::vector<double> result = U.ExtractColumn(m_nSteps);
    std::transform(result.begin(), result.end(), result.begin(), [discount](double x) { return discount * x; });
    return PricerOutput(result, 0.99);
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

static void testGeneratePath()
{
    double date = 1.0;
    BSMModel model(100, 0.1, 0.05, 0, date); // Low volatility to have narrow final prices sample

    std::size_t nPaths = 10000;
    std::size_t nSteps = 5;
    MonteCarloPricer mc(nPaths, nSteps);
    double maturity = 2.0;

    Matrix pricePaths = mc.generatePricePaths(model, maturity);
    std::vector<double> finalPrices = pricePaths.ExtractColumn(5);
    ASSERT_APPROX_EQUAL(mean(finalPrices), exp(model.InterestRate() * (maturity - model.Date())) * model.StockPrice(), 0.5);
    
    MonteCarloPricer mc2(1, nSteps); // 1 path, 10 000 steps
    std::vector<double> paths = (mc2.generatePricePaths(model, maturity)).asRow();
    std::vector<double> dates = linspace(0, maturity, nSteps + 1);
    LineChart lineChart("Stock price path");
    lineChart.setSeries(dates, paths);
    lineChart.writeAsHTML("examplePricePath.html");
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

/*
 Numerical example from Paul Wilmott's book
 */
/*
static void testLSWeights()
{
    Matrix pricePaths("100,92.30759,107.7357,98.04343,110.7416,95.34586;100,103.4446,111.9465,110.9322,117.8379,119.4419;100,111.2298,121.2417,137.1683,145.1687,133.1789;100,105.7152,115.0572,99.73054,100.6804,100.9471;100,98.47278,96.5825,91.32007,80.63689,82.1163;100,94.40168,94.16078,87.83702,93.84797,93.45847;100,106.7042,125.264,139.4822,132.0177,126.2041;100,84.37568,76.60055,76.21345,80.85454,95.19434;100,94.21698,88.00477,90.81541,88.63676,84.80556;100,99.81029,105.2631,101.747,103.1483,107.3703");
    BSMModel model(100, 0.2, 0.05);
    std::shared_ptr<VanillaOption> USPut = VanillaOption::newOption(100, 1, put, American);
    PolynomialBasis basis(2);
    MonteCarloPricer pricer(10, 5, &basis);
    Matrix result = pricer.LSWeights(*USPut, model, pricePaths);
    Matrix expected("903.58371,579.009566,-203.94825,-1075.2312;-18.505142,-12.62883,5.66129963,25.3466805;0.0955491,0.07027471,-0.0360753,-0.1472459");
    result.assertEquals(expected, 1e-4);
}
 */
// As test passed, we will make this method private

static void testUSOptionPricing()
{
    LaguerreBasis basis(2);
    MonteCarloPricer mc(10, 5, &basis);
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
    //TEST(testLSWeights);
    TEST(testUSOptionPricing);
}
