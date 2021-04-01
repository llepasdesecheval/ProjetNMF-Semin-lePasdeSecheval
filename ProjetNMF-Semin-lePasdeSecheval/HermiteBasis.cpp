//
//  HermiteBasis.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "HermiteBasis.hpp"

HermiteBasis::HermiteBasis(std::size_t nbRegressors) : Basis(nbRegressors) {}

double HermiteBasis::getKthRegressor(double x, std::size_t k) const
{
    if (k == 0)
    {
        return 1;
    }
    else
    {
        double result = 0;
        double parity = 1;
        double power = std::pow(x, k);
        for (unsigned int i = 0; i < k/2 + 1; ++i)
        {
            double a = static_cast<double>(factorial(k))/(factorial(i) * factorial(k - 2 * i)) * parity * power;
            result += a;
            parity *= - 0.5;
            power /= x * x;
        }
        return result;
    }
}

static void testHermiteEvaluation()
{
    double x = 3.15;
    HermiteBasis test(6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 0), 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 1), x, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 2), x * x - 1, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 3), (x * x - 3) * x, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 4), (x * x - 6) * x * x + 3, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 5), ((x * x - 10) * x * x + 15) * x, 1e-6);
    ASSERT_APPROX_EQUAL(test.getKthRegressor(x, 6), ((x * x - 15) * x * x + 45) * x * x - 15, 1e-6);
}

void testHermiteBasis()
{
    TEST(testHermiteEvaluation);
}
