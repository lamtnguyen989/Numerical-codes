#ifndef INTEGRATION_DATA_H
#define INTEGRATION_DATA_H

#include <array>
#include <functional>

class Integration_data
{
    public:
        // Constructor signatures
        Integration_data(std::function<double(double)> f, std::array<double, 2> range, double result);
        Integration_data(std::function<double(double)> f, double upper_bound, double lower_bound, double result);

        // Getters
        std::function<double(double)> get_integrand()   {return integrand;};
        std::array<double, 2> get_integration_bounds()  {return integration_range;}
        double get_analytical_result()                  {return analytical_result;};
        double get_lower_integration_bound()            {return integration_range[0];}
        double get_upper_integration_bound()            {return integration_range[1];}
    private:
        std::function<double(double)>   integrand;
        std::array<double, 2>           integration_range;
        double                          analytical_result;
};

Integration_data::Integration_data(std::function<double(double)> f, std::array<double, 2> range, double result)
    : integrand(f)
    , integration_range(range)
    , analytical_result(result)
{}

Integration_data::Integration_data(std::function<double(double)> f, double lower_bound, double upper_bound, double result)
    : integrand(f)
    , analytical_result(result)
{
    integration_range = {lower_bound, upper_bound};
}

#endif // INTEGRATION_DATA_H