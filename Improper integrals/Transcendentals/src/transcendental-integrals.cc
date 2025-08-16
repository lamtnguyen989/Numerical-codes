#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "integration_data.h"   // For the object that encapsulate integration data

using namespace boost::math::quadrature;

int main()
{
    // Initializing integrands
    std::function<double(double)> f1 = [](double x) {return (std::cos(x) - 1) / std::pow(x,2);};    // (cos(x) - 1) / x^2
    std::function<double(double)> f2 = [](double x) {return std::log(x) / (std::pow(x,4) + 1);};    // ln(x) / (x^4-1)

    // Finalizing integration datum
    std::vector<Integration_data> datum;
    datum.push_back(Integration_data(f1, {0, std::numeric_limits<double>::infinity()}, -M_PI/2));
    datum.push_back(Integration_data(f2, 0, std::numeric_limits<double>::infinity(), -(std::pow(M_PI,2)*std::sqrt(2))/16));

    // Integration loop
    double estimated_error;
    double result;
    for (unsigned int i = 0; i < datum.size(); i++)
    {
        result = gauss_kronrod<double, 61>::integrate(datum[i].get_integrand(), datum[i].get_lower_integration_bound(), datum[i].get_upper_integration_bound(), 12, 1e-14, &estimated_error);
        std::cout << "Actual integration result: " << std::setprecision(10) << result << std::endl;
        std::cout << "Estimated error: " << estimated_error << std::endl;
        if (!std::isnan(datum[i].get_analytical_result()))
        {
            std::cout << "Actual error: " << std::abs(datum[i].get_analytical_result() - result) << std::endl << std::endl;
        }
    }
}