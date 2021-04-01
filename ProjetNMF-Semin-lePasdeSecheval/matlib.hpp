//
//  matlib.hpp
//  ProjetNMF-Semin-lePasdeSecheval
//
//  Created by Louis le Pas de Sécheval on 01/04/2021.
//

#if !defined MATLIB_HPP
#define MATLIB_HPP

#include "stdafx.hpp"
#include "Basis.hpp"

/** Pi */
const double pi = 3.1415926535897932384626433832795028842;

/**
 Creates a linearly spaced vector.
 @param from Starting point
 @param to Stopping point, included in the output
 @param numPoints Number of points, including start and stop
 */
std::vector<double> linspace(double from, double to, std::size_t numPoints);

/////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Standard Gaussian distribution functions
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 Standard Gaussian distribution.
 @param x Computation point
 @param cumulative If true, computes the cumulative distribution function; else, computes the probability density function.
 */
double NormSDist(double x, bool cumulative = true);

/**
 Computes the inverse cumulative distribution function of the standard Gaussian distribution using Moro's algorithm
 */
double NormSDistInv(double x);

/////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Statistical functions
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 Finds the minimum of a vector.
 @param data Input vector
 */
double min(const std::vector<double>& data);

/**
 Finds the maximum of a vector.
 @param data Input vector
 */
double max(const std::vector<double>& data);

/**
 Computes the sum of a vector
 @param data Input data
 */
double sum(const std::vector<double>& data);

/**
Computes the mean of a vector
@param data Input data
*/
double mean(const std::vector<double>& data);

/**
 Computes the variance of a vector.
 @param data Input data
 @param sample If true, returns the sample estimate of variance; else, return the population variance
*/
double variance(const std::vector<double>& data, bool sample = true);

/**
 Computes the standard deviation of a vector.
 @param data Input data
 @param sample If true, returns the sample estimate of standard deviation; else, return the population standard deviation
*/
double stDev(const std::vector<double>& data, bool sample = true);

/**
 Computes the covariance of two vectors.
 @param x First series of input data
 @param y Second series of input data
 @param sample If true, returns the sample estimate of covariance; else, return the population covariance
*/
double covariance(const std::vector<double>& x, const std::vector<double>& y, bool sample = true);

/**
 Computes the correlation of two vectors.
 @param x First series of input data
 @param y Second series of input data
*/
double correlation(const std::vector<double>& x, const std::vector<double>& y);

/////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Random number simulations
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/**
Creates uniformly distributed random numbers in a range
@param n Size of output
@param a Lower bound
@param b Upper bound
*/
std::vector<double> randUniform(std::size_t n, double a = 0, double b = 1);

/**
 Creates Gaussian random numbers using the Marsaglia polar method
 @param n Size of output
 @param sigma Standard deviation
 @param mu Mean of the distribution
 */
std::vector<double> randNMarsaglia(std::size_t n, double sigma = 1, double mu = 0);

/**
 Creates Gaussian random numbers using the Box-Muller method
 @param n Size of output
 @param sigma Standard deviation
 @param mu Mean of the distribution
 */
std::vector<double> randNBoxMuller(std::size_t n, double sigma = 1, double mu = 0);

/**
Creates Gaussian random numbers using the inverse CDF
@param n Size of output
@param sigma Standard deviation
@param mu Mean of the distribution
*/
std::vector<double> randNInversion(std::size_t n, double sigma = 1, double mu = 0);

/**
Creates Gaussian random numbers using standard implementation
@param n Size of output
@param sigma Standard deviation
@param mu Mean of the distribution
*/
std::vector<double> randN(std::size_t n, double sigma = 1, double mu = 0);

/** Seeds the default random number generator */
void setDefaultSeed();

/////////////////////////////////////////////////////////////////////////////////////////////////
///
/// Regression
///
/////////////////////////////////////////////////////////////////////////////////////////////////

/**
 Returns the coefficients of a regression for a given set of basis functions
 @param X Vector of input
 @param Y Vector of output
 @param basis Basis functions to use
 */
std::vector<double> RegressionCoefficients(const std::vector<double>& X, const std::vector<double>& Y,
                                           const Basis& basis);

unsigned long long factorial(unsigned long long n);

/** Test function for matlib.hpp */
void testMatlib();

#endif /* matlib_hpp */
