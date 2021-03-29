//
//  PricerOutput.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 29/03/2021.
//

#if !defined PRICEROUTPUT_HPP
#define PRICEROUTPUT_HPP

#include "stdafx.hpp"

class PricerOutput
{
public:
    /**
    Constructor from a Monte Carlo simulation sample
    @param sample 1-dimensional data vector
    @param confidenceLevel Desired confidence level, between 0 and 1.
    */
    PricerOutput(std::vector<double> sample, double confidenceLevel);

    /**
    Constructor from a known estimate and confidence interval size
    @param estimate Estimate
    @param confidenceRadius Confidence interval radius
    */
    PricerOutput(double estimate, double confidenceRadius);

    /**
    Constructor from a closed-form formula; confidence interval has width 0.
    @param estimate Estimate
    */
    PricerOutput(double estimate);

    double Estimate();
    double Confidence();

private:
    double m_estimate;
    double m_confidence;
};

#endif /* PricerOutput_hpp */
