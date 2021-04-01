//
//  matlib.cpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#include "matlib.hpp"
#include "PolynomialBasis.hpp"
#include <limits>

const double invroot2pi = 1.0/std::sqrt(2 * pi);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Horner polynomials (for Gaussian inverse CDF)
///
/// Horner's algorithm allows faster computation of polynomials
/// Complexity O(n) instead of O(n^2), n highest polynomial order
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static inline double hornerFunction(double x, double a0, double a1)
{
    return a0 + x * a1;
}

static inline double hornerFunction(double x, double a0, double a1, double a2)
{
    return a0 + x * hornerFunction(x, a1, a2);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3)
{
    return a0 + x * hornerFunction(x, a1, a2, a3);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7);
}

static inline double hornerFunction(double x, double a0, double a1, double a2, double a3, double a4, double a5, double a6, double a7, double a8)
{
    return a0 + x * hornerFunction(x, a1, a2, a3, a4, a5, a6, a7, a8);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Constants required for Moro's algorithm
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static const double a0 = 2.50662823884;
static const double a1 = -18.61500062529;
static const double a2 = 41.39119773534;
static const double a3 = -25.44106049637;
static const double b1 = -8.47351093090;
static const double b2 = 23.08336743743;
static const double b3 = -21.06224101826;
static const double b4 = 3.13082909833;
static const double c0 = 0.3374754822726147;
static const double c1 = 0.9761690190917186;
static const double c2 = 0.1607979714918209;
static const double c3 = 0.0276438810333863;
static const double c4 = 0.0038405729373609;
static const double c5 = 0.0003951896511919;
static const double c6 = 0.0000321767881768;
static const double c7 = 0.0000002888167364;
static const double c8 = 0.0000003960315187;

std::vector<double> linspace(double from, double to, std::size_t numPoints)
{
    ASSERT(numPoints >= 2);
    std::vector<double> ret(numPoints, 0.0);
    double step = (to - from)/(numPoints - 1);
    double current = from;
    for (std::size_t i = 0; i < numPoints; ++i)
    {
        ret[i] = current;
        current += step;
    }
    return ret;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Standard Gaussian distribution functions
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double NormSDist(double x, bool cumulative)
{
    if (cumulative) {
        DEBUG_PRINT("normCDF(" << x << ")");
        return std::erfc(-x / sqrt(2)) / 2;
    }
    else {
        DEBUG_PRINT("normPDF(" << x << ")");
        return invroot2pi * exp(-0.5 * x * x);
    }
}

double NormSDistInv(double x)
{
    if ((x < 0.0) || (x > 1.0))
    {
        throw std::logic_error("Input must be between 0.0 and 1.0");
    }
    if (x == 0)
        return -std::numeric_limits<double>::infinity();
    else if (x == 1)
        return std::numeric_limits<double>::infinity();
    else
    {
        DEBUG_PRINT("normInv(" << x << ")");
        double y = x - 0.5;
        double r;
        double s;
        double t;

        if (abs(y) < 0.42)
        {
            r = y * y;
            DEBUG_PRINT("Case 1, r=" << r);
            return y * hornerFunction(r, a0, a1, a2, a3) / hornerFunction(r, 1.0, b1, b2, b3, b4);
        }
        else
            r = (y < 0.0) ? x : 1.0 - x;

        DEBUG_PRINT("Case 2, r=" << r);
        s = log(-log(r));
        t = hornerFunction(s, c0, c1, c2, c3, c4, c5, c6, c7, c8);
        return (x > 0.5) ? t : -t;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Statistical functions
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double min(const std::vector<double>& data)
{
    ASSERT(data.size() > 0);
    return *std::min_element(data.begin(), data.end());
}

double max(const std::vector<double>& data)
{
    ASSERT(data.size() > 0);
    return *std::max_element(data.begin(), data.end());
}

double sum(const std::vector<double>& data)
{
    ASSERT(data.size() > 0);
    return std::accumulate(data.begin(), data.end(), 0.);
}

double mean(const std::vector<double>& data)
{
    // The test for size is already in sum;
    return sum(data) / data.size();
}

double variance(const std::vector<double>& data, bool sample)
{
    std::size_t n = data.size();
    ASSERT(n > 0);
    double sum = 0;
    double sumSquared = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double elt = data[i];
        sum += elt;
        sumSquared += elt * elt;
    }
    /*
    Why all this? Performance.
    Had we called the mean function inside our variance function,
    we would have implicitly looped twice!

    Here, we only loop once.
    */
    return (sumSquared - sum * sum / n) / (sample ? n - 1 : n);
}

double stDev(const std::vector<double>& data, bool sample)
{
    return sqrt(variance(data, sample));
}

double covariance(const std::vector<double>& x, const std::vector<double>& y, bool sample)
{
    std::size_t n = x.size();
    ASSERT(n > 0);
    ASSERT(y.size() == n);
    double sum_x = 0;
    double sum_y = 0;
    double sumCross = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double elt_x = x[i];
        double elt_y = y[i];
        sum_x += elt_x;
        sum_y += elt_y;
        sumCross += elt_x * elt_y;
    }
    return (sumCross - sum_x * sum_y / n) / (sample ? n - 1 : n);
}

double correlation(const std::vector<double>& x, const std::vector<double>& y)
{
    std::size_t n = x.size();
    ASSERT(n > 0);
    ASSERT(y.size() == n);
    double sum_x = 0;
    double sum_y = 0;
    double sumSq_x = 0;
    double sumSq_y = 0;
    double sumCross = 0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double elt_x = x[i];
        double elt_y = y[i];
        sum_x += elt_x;
        sum_y += elt_y;
        sumSq_x += elt_x * elt_x;
        sumSq_y += elt_y * elt_y;
        sumCross += elt_x * elt_y;
    }
    return (sumCross - sum_x * sum_y / n) / sqrt((sumSq_x - sum_x * sum_x / n) * (sumSq_y - sum_y * sum_y / n));
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Random Number Simulation
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static std::random_device rand_dev;
static std::mt19937_64 mersenneTwister(rand_dev());

std::vector<double> randUniform(std::size_t n, double a, double b)
{
    std::uniform_real_distribution<> dis(a, b);
    std::vector<double> result(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        result[i] = dis(mersenneTwister);
    }
    return result;
}

std::vector<double> randNMarsaglia(std::size_t n, double sigma, double mu)
{
    ASSERT(sigma >= 0);
    std::vector<double> result(n);
    if (sigma == 0)
    {
        std::fill(result.begin(), result.end(), mu);
    }
    else
    {
        std::size_t roundUpN = n;
        if (n % 2 != 0)
        {
            roundUpN++;
            result.push_back(0);
        }
        result = randUniform(roundUpN, -1.0, 1.0);
        for (std::size_t i = 0; i < roundUpN; i += 2)
        {
            double u1 = result[i];
            double u2 = result[i + 1];
            double s = u1 * u1 + u2 * u2;
            while ((s >= 1) || (s <= 0))
            {
                std::vector<double> tryAgain = randUniform(2, -1.0, 1.0);
                u1 = tryAgain[0];
                u2 = tryAgain[1];
                s = u1 * u1 + u2 * u2;
            }
            double r = std::sqrt(-2 * log(s)/s);
            double n1 = mu + sigma * r * u1;
            double n2 = mu + sigma * r * u2;
            result[i] = n1;
            result[i + 1] = n2;
        }
        if (n != roundUpN)
        {
            result.pop_back(); // if we have one more element than requested remove it
        }
    }
    return result;
}

std::vector<double> randNBoxMuller(std::size_t n, double sigma, double mu)
{
    ASSERT(sigma >= 0);
    std::vector<double> result(n);
    if (sigma == 0)
    {
        std::fill(result.begin(), result.end(), mu);
    }
    else
    {
        std::size_t roundUpN = n;
        if (n % 2!= 0)
        {
            roundUpN++;
            result.push_back(0);
        }
        result = randUniform(roundUpN);
        for (std::size_t i = 0; i < roundUpN; i += 2)
        {
            double u1 = result[i];
            double u2 = result[i+1];
            double r = std::sqrt(-2 * log(u1));
            double theta = 2 * pi * u2;
            double n1 = mu + sigma * r * cos(theta);
            double n2 = mu + sigma * r * sin(theta);
            result[i] = n1;
            result[i + 1] = n2;
        }
        if (n != roundUpN)
        {
            result.pop_back(); // if we have one more element than requested remove it
        }
    }
    return result;
}

std::vector<double> randNInversion(std::size_t n, double sigma, double mu)
{
    ASSERT(sigma >= 0);
    std::vector<double> result(n);
    if (sigma == 0)
    {
        std::fill(result.begin(), result.end(), mu);
    }
    else
    {
        result = randUniform(n);
        auto GaussInv = [mu, sigma](double x) {return mu + sigma * NormSDistInv(x); };
        std::transform(result.begin(), result.end(), result.begin(), GaussInv);
    }
    return result;
}

std::vector<double> randN(std::size_t n, double sigma, double mu)
{
    // This function got the "simple" name because it is the fastest one, hence we will use it.
    ASSERT(sigma >= 0);
    std::size_t roundUpN = (n % 2 == 0) ? n : n + 1;
    std::vector<double> result(roundUpN);
    if (sigma == 0)
    {
        std::fill(result.begin(), result.end(), mu);
    }
    else
    {
        std::normal_distribution<> dis(mu, sigma);
        for (std::size_t i = 0; i < roundUpN; i += 2)
        {
            double sample = dis(mersenneTwister);
            result[i] = sample;
            result[i + 1] = 2 * mu - sample;
        }
    }
    while (result.size() > n)
    {
        result.pop_back();
    }
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Regression
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> RegressionCoefficients(const std::vector<double>& X, const std::vector<double>& Y,
                                           const Basis& basis)
{
    Matrix L = basis.RegressorMatrix(X);
    /*
     Our regression writes Y = L * beta + error, where L is the regressor matrix.
     The OLS solution to this is beta = (L'L)^-1 L' Y.
     */
    Matrix A = L.transpose() * L;
    std::vector<double> beta = A.Inverse() * L.transpose() * Y;
    return beta;
}

unsigned long long factorial(unsigned long long n)
{
    return (n > 1) ? n * factorial(n - 1) : 1;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Tests
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void testNormSDist()
{
    // Test bounds
    ASSERT(NormSDist(0.3) > 0);
    ASSERT(NormSDist(0.3) < 1);
    // Test extreme values
    ASSERT_APPROX_EQUAL(NormSDist(-1e10), 0, 0.001);
    ASSERT_APPROX_EQUAL(NormSDist(1e10), 1.0, 0.001);
    // Test increasing
    ASSERT(NormSDist(0.3) < NormSDist(0.5));
    // Test symmetry
    ASSERT_APPROX_EQUAL(NormSDist(0.3), 1 - NormSDist(-0.3), 0.0001);
    ASSERT_APPROX_EQUAL(NormSDist(0.0), 0.5, 0.0001);
    // Test inverse
    ASSERT_APPROX_EQUAL(NormSDist(NormSDistInv(0.3)), 0.3, 0.0001);
    // Test well-known value
    ASSERT_APPROX_EQUAL(NormSDist(1.96), 0.975, 0.001);
}

static void testNormInv()
{
    ASSERT_APPROX_EQUAL(NormSDistInv(0.975), 1.96, 0.01);
}

static std::vector<double> testVector = { 1, 5, 3, 9, 7 };
static std::vector<double> testVector2 = { 2, 4, 0, 6, 8 };

static void testMin() {
    ASSERT_APPROX_EQUAL(min(testVector), 1.0, 0.001);
    ASSERT_APPROX_EQUAL(min(testVector2), 0.0, 0.001);
}

static void testMax() {
    ASSERT_APPROX_EQUAL(max(testVector), 9.0, 0.001);
    ASSERT_APPROX_EQUAL(max(testVector2), 8.0, 0.001);
}

void testMean()
{
    ASSERT_APPROX_EQUAL(mean(testVector), 5.0, 0.0001);
}

void testStDev()
{
    ASSERT_APPROX_EQUAL(stDev(testVector), 3.1623, 0.0001);
    ASSERT_APPROX_EQUAL(stDev(testVector, false), 2.8284, 0.0001);
}

void testCovariance()
{
    ASSERT_APPROX_EQUAL(covariance(testVector, testVector2), 8, 0.0001);
    ASSERT_APPROX_EQUAL(covariance(testVector, testVector2, false), 6.4, 0.0001);
    // Covariance(X, X) = Variance(X)
    ASSERT_APPROX_EQUAL(covariance(testVector, testVector), variance(testVector), 0.001);
    ASSERT_APPROX_EQUAL(covariance(testVector, testVector, false), variance(testVector, false), 0.001);
}

void testCorrelation()
{
    ASSERT_APPROX_EQUAL(correlation(testVector, testVector2), 0.8, 0.0001);
    // Correlation(X, X) = 1
    ASSERT_APPROX_EQUAL(correlation(testVector, testVector), 1.0, 0.0001);
}

void setDefaultSeed()
{
    mersenneTwister.seed(std::mt19937_64::default_seed);
}

static void testRandUniform()
{
    setDefaultSeed();
    std::vector<double> v = randUniform(100000, -5, 4);
    ASSERT(v.size() == 100000);
    ASSERT_APPROX_EQUAL(mean(v), -0.5, 0.01);
    ASSERT(max(v) < 4.0);
    ASSERT(min(v) > -5.0);
}

static void testBoxMuller()
{
    setDefaultSeed();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> eps = randNBoxMuller(1000000);
    auto t2 = std::chrono::high_resolution_clock::now();
    ASSERT_APPROX_EQUAL(mean(eps), 0.0, 0.01);
    ASSERT_APPROX_EQUAL(stDev(eps), 1.0, 0.01);
    
    // sizing code is the hardest so let's check it
    std::vector<double> v1 = randNBoxMuller(3);
    ASSERT(v1.size() == 3);
    std::vector<double> v2 = randNBoxMuller(4);
    ASSERT(v2.size() == 4);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Number of simulations: " << eps.size() << ".\n";
    std::cout << "Execution time: " << elapsed << " ms." << std::endl;
}

static void testMarsaglia()
{
    setDefaultSeed();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> eps = randNMarsaglia(1000000);
    auto t2 = std::chrono::high_resolution_clock::now();
    ASSERT_APPROX_EQUAL(mean(eps), 0.0, 0.01);
    ASSERT_APPROX_EQUAL(stDev(eps), 1.0, 0.01);
    
    // sizing code is the hardest so let's check it
    std::vector<double> v1 = randNMarsaglia(3);
    ASSERT(v1.size() == 3);
    std::vector<double> v2 = randNMarsaglia(4);
    ASSERT(v2.size() == 4);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Number of simulations: " << eps.size() << ".\n";
    std::cout << "Execution time: " << elapsed << " ms." << std::endl;
}

static void testGaussInversion() {
    setDefaultSeed();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> eps = randNInversion(1000000);
    auto t2 = std::chrono::high_resolution_clock::now();
    ASSERT_APPROX_EQUAL(mean(eps), 0.0, 0.01);
    ASSERT_APPROX_EQUAL(stDev(eps), 1.0, 0.01);
    
    // sizing code is the hardest so let's check it
    std::vector<double> v1 = randNInversion(3);
    ASSERT(v1.size() == 3);
    std::vector<double> v2 = randNInversion(4);
    ASSERT(v2.size() == 4);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Number of simulations: " << eps.size() << ".\n";
    std::cout << "Execution time: " << elapsed << " ms." << std::endl;
}

static void testGaussStandard() {
    setDefaultSeed();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::vector<double> eps = randN(1000000);
    auto t2 = std::chrono::high_resolution_clock::now();
    ASSERT_APPROX_EQUAL(mean(eps), 0.0, 0.01);
    ASSERT_APPROX_EQUAL(stDev(eps), 1.0, 0.01);
    
    eps = randN(1000000, 4.2, 3.0);
    ASSERT_APPROX_EQUAL(mean(eps), 3.0, 0.01);
    ASSERT_APPROX_EQUAL(stDev(eps), 4.2, 0.01);
    
    // sizing code is the hardest so let's check it
    std::vector<double> v1 = randN(3);
    ASSERT(v1.size() == 3);
    std::vector<double> v2 = randN(4);
    ASSERT(v2.size() == 4);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    std::cout << "Number of simulations: " << eps.size() << ".\n";
    std::cout << "Execution time: " << elapsed << " ms." << std::endl;
}

static void testPolynomialRegression()
{
    /*
     We test a polynomial regression of degree 2 carried in Excel
     The generating equation was Y = X^2 + 3X + eps, where eps was a standard Gaussian white noise
     */
    std::vector<double> X = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                            21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
    std::vector<double> Y = {0.213570929746912, 3.36681662053155, 8.51791732086752, 15.9663982724039, 27.7592776859739,
                            38.233723126632, 53.0996525516416, 69.5968081484078, 87.3637826102275, 108.277268247296,
                            130.68154619144, 153.13353314145, 179.150012621737, 206.697431675618, 237.652416222334,
                            270.8072491269, 304.915317879837, 341.892768212285, 377.469531547685, 417.774267092854,
                            460.011118384517, 503.28089352627, 548.520995953777, 596.192471619217, 647.487317547541,
                            699.749590077496, 753.108262939826, 810.032140426296, 867.456891346427, 927.289250950271,
                            991.435137390862};
    std::vector<double> coefSq = {-0.919333, 3.061561, 0.998626}; // Coefficients yielded by Excel
    PolynomialBasis test(2);
    std::vector<double> regressionSqOutput = RegressionCoefficients(X, Y, test);
    ASSERT(regressionSqOutput.size() == coefSq.size());
    for (std::size_t i = 0; i < coefSq.size(); ++i)
    {
        ASSERT_APPROX_EQUAL(coefSq[i], regressionSqOutput[i], 1e-6);
    }
    
    /**
     Let us try a linear regression on the same dataset
     */
    std::vector<double> coefLin = {-145.720074, 33.020335};
    PolynomialBasis testLin(1);
    std::vector<double> regressionLinOutput = RegressionCoefficients(X, Y, testLin);
    ASSERT(regressionLinOutput.size() == coefLin.size());
    for (std::size_t i = 0; i < coefLin.size(); ++i)
    {
        ASSERT_APPROX_EQUAL(coefLin[i], regressionLinOutput[i], 1e-6);
    }
}

static void testFactorial()
{
    std::vector<int> v = {0, 1, 2, 3, 4, 5, 6};
    std::transform(v.begin(), v.end(), v.begin(), &factorial);
    std::vector<int> expected = {1, 1, 2, 6, 24, 120, 720};
    ASSERT(v == expected);
}

#pragma endregion

void testMatlib() {
    TEST(testNormSDist);
    TEST(testNormInv);
    TEST(testMin);
    TEST(testMax);
    TEST(testMean);
    TEST(testStDev);
    TEST(testCovariance);
    TEST(testCorrelation);
    TEST(testRandUniform);
    TEST(testMarsaglia);
    TEST(testBoxMuller);
    TEST(testGaussInversion);
    TEST(testGaussStandard);
    TEST(testPolynomialRegression);
    TEST(testFactorial);
}
