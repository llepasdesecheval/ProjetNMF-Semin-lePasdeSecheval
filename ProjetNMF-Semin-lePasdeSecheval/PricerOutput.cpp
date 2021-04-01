//
//  PricerOutput.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 29/03/2021.
//

#include "PricerOutput.hpp"
#include "matlib.hpp"

PricerOutput::PricerOutput(std::vector<double> sample, double confidenceLevel) : m_estimate(mean(sample))
{
    ASSERT(sample.size() > 0);
    ASSERT((confidenceLevel > 0) && (confidenceLevel < 1));
    m_confidence = NormSDistInv((1 + confidenceLevel) / 2) * stDev(sample) / sqrt(sample.size());
}

PricerOutput::PricerOutput(double estimate, double confidenceRadius) : m_estimate(estimate), m_confidence(confidenceRadius) {}

PricerOutput::PricerOutput(double estimate) : m_estimate(estimate), m_confidence(0) {}

double PricerOutput::Estimate()
{
    return m_estimate;
}

double PricerOutput::Confidence()
{
    return m_confidence;
}
