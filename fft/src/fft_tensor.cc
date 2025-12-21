#include <iostream>
#include <complex>
#include <vector>

#include "tensor.h"

using Complex = std::complex<float>;

int main()
{
    Tensor<Complex> tensor({1.0}, {1});
    std::cout << "Hello world " << tensor.order() <<std::endl;
}