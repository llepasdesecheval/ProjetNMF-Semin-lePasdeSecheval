//
//  MonteCarloPricer.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined MONTECARLOPRICER_HPP
#define MONTECARLOPRICER_HPP

#include "Pricer.hpp"
#include "Basis.hpp"
#include "Matrix.hpp"

class MonteCarloPricer : public Pricer
{
public:
    /**
     Constructor
     @param nScenarii Number of simulated paths
     @param nSteps Number of steps in each path
     */
    MonteCarloPricer(std::size_t nScenarii, std::size_t nSteps, Basis* basis = nullptr);

    /**
     Stock price path generator. Rows are scenarios, columns are time steps.
     @param model Black-Scholes-Merton model (BSMModel) containing the simulation parameters
     @param toDate Time horizon of the simulations
     */
    Matrix generatePricePaths(const BSMModel& model, double toDate) const;
    
    /**
     For Monte Carlo pricing, we chose to estimate the price with a 99% confidence level
     */
    PricerOutput price(const VanillaOption& option, const BSMModel& model) const override;
    
private:
    std::size_t m_nScenarii;
    std::size_t m_nSteps;
    Basis* m_basis;
    
    /**
     Nested Monte Carlo simulation (1-step, hence returns a vector)
     @param St Spot price
     @param horizon Nested simulations horizon
     @param nSuccessors Number of nested simulations
     @param model Black-Scholes-Merton model (BSMModel) containing the simulation parameters
     */
    std::vector<double> nestedMC(double St, double horizon, std::size_t nSuccessors, const BSMModel& model) const;
    
    PricerOutput LSPrice(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const;
    PricerOutput ABPrice(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const;
    Matrix LSWeights(const VanillaOption& option, const BSMModel& model, const Matrix& pricePaths) const;
};

void testMonteCarloPricer();

#endif /* MonteCarloPricer_hpp */
