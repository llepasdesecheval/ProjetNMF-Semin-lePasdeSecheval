//
//  LaguerreBasis.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "LaguerreBasis.hpp"

LaguerreBasis::LaguerreBasis(std::size_t nbRegressors) : Basis(nbRegressors) {}

double LaguerreBasis::getKthRegressor(double x, std::size_t k) const
{
    double result = 0;
    double power = 1;
    for (unsigned int i = 0; i < k + 1; ++i)
    {
        result += static_cast<double>(factorial(k))/(factorial(i) * factorial(i) * factorial(k - i)) * power;
        power *= -x;
    }
    return result;
}

static void testLaguerreEvaluation()
{
    double x = 3.15;
    LaguerreBasis test(6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 0), 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 1), 1 - x, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 2), (0.5 * x - 2) * x + 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 3), ((- 1.0/6 * x + 1.5) * x - 3) * x + 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 4), (((1.0/24 * x - 2.0/3) * x + 3) * x - 4) * x + 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 5), ((((- 1.0/120 * x + 5.0/24) * x - 5.0/3) * x + 5) * x - 5) * x + 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 6), (((((1.0/720 * x - 0.05) * x + 0.625) * x - 10.0/3) * x + 7.5) * x - 6) * x + 1, 1e-6);
}

void testLaguerreBasis()
{
    TEST(testLaguerreEvaluation);
}
