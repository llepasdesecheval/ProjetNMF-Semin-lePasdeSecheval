//
//  PolynomialBasis.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 29/03/2021.
//

#include "PolynomialBasis.hpp"

PolynomialBasis::PolynomialBasis(std::size_t nbRegressors) : Basis(nbRegressors) {}

double PolynomialBasis::getKthRegressor(double x, std::size_t k) const
{
    ASSERT(k <= m_nbRegressors);
    double result = 1;
    for (int i = 0; i < k; ++i)
    {
        result *= x;
    }
    return result;
}

static void testPolynomialEvaluation()
{
    PolynomialBasis test(1);
    double x = 1345468.8976525;
    ASSERT(test.getKthRegressor(x, 0) == 1);
    ASSERT(test.getKthRegressor(x, 1) == x);
}

static void testRegressorMatrix()
{
    std::vector<double> testVector{ 1.0, 2.0, 3.0, 4.0, 5.0 };
    std::size_t n = testVector.size();
    PolynomialBasis test(4);
    Matrix output = test.RegressorMatrix(testVector);
    ASSERT(output.nCols() == 5);
    ASSERT(output.nRows() == n);
    printMatrix(output);
    Matrix expected("1,1,1,1,1;1,2,4,8,16;1,3,9,27,81;1,4,16,64,256;1,5,25,125,625");
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < 5; ++j)
        {
            ASSERT_APPROX_EQUAL(output(i, j), expected(i, j), 1e-6);
        }
    }
}

static void testRegressorVector()
{
    double x = 8.735957312;
    PolynomialBasis test(4);
    std::vector<double> output = test.RegressorVector(x);
    std::vector<double> expected { 1.0, 8.735957312, 76.3169501570863, 666.701618754337, 5824.27688127919 };
    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_APPROX_EQUAL(output[i], expected[i], 1e-6);
    }
}

void testPolynomialBasis()
{
    TEST(testPolynomialEvaluation);
    TEST(testRegressorMatrix);
    TEST(testRegressorVector);
}
