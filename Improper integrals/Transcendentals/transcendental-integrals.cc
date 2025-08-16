#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using namespace boost::math::quadrature;

// integrand: (cos(x) - 1) / x^2
auto f1 = [](double x) {return (std::cos(x) - 1) / std::pow(x,2);};

// integrand: ln(x) / (x^4-1)
auto f2 = [](double x) {return std::log(x) / (std::pow(x,4) + 1);};


// Containers for integrands and analytical results
std::vector<std::function<double(double)>> integrands = {f1, f2};
std::vector<double> analytical_results = {-M_PI/2, -(std::pow(M_PI,2)*std::sqrt(2))/16};

int main()
{
    double estimated_error;
    double result;

    for (unsigned int i = 0; i < integrands.size(); i++)
    {
        result = gauss_kronrod<double, 61>::integrate(integrands[i], 0, std::numeric_limits<double>::infinity(), 12, 1e-14, &estimated_error);
        std::cout << "Actual integration result: " << std::setprecision(10) << result << std::endl;
        std::cout << "Estimated error: " << estimated_error << std::endl;
        std::cout << "Actual error: " << std::abs(analytical_results[i] - result) << std::endl << std::endl;
    }
}