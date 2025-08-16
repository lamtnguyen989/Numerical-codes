#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace boost::math::quadrature;

// integrand: e^{-x^2}
auto standard_gaussian = [](double x) {return std::exp(-std::pow(x,2));};

int main()
{
    // Estimated error container
    double estimated_error;

    // Computing the integral
    double result = gauss_kronrod<double, 61>::integrate(standard_gaussian, -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), 6, 1e-14, &estimated_error);

    // Print
    std::cout << "Integration of standard Gaussian over the whole real line" << std::endl << std::endl;
    std::cout << "Actual integration result: " << std::setprecision(result) << std::endl;
    std::cout << "Estimated error: " << estimated_error << std::endl;
    std::cout << "Actual error: " << std::abs(std::sqrt(M_PI) - result) << std::endl;
}