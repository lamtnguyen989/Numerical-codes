#include <iostream>
#include <complex>
#include <vector>

void fft(std::vector<std::complex<double>> &X)
{
    // Length and helper vectors
    const unsigned int N = X.size();
    std::vector<std::complex<double>> evens(N/2);
    std::vector<std::complex<double>> odds(N/2);

    if (N == 1 )
        return;     // Recursion base case

    // Copy the evens and odds
    for (unsigned int k = 0; k < N/2; k++) {
        evens.at(k) = X.at(2*k);
        odds.at(k) = X.at(2*k + 1);
    }

    // Recurse the DFT
    fft(evens);
    fft(odds);

    // Computation result
    for (unsigned int k = 0; k < N/2; k++) {
        std::complex<double> root_of_unity = std::polar(1.0, -2*M_PI*k / N);
        X.at(k) = evens.at(k) + root_of_unity * odds.at(k);
        X.at(k + N/2) = evens.at(k) - root_of_unity * odds.at(k);
    }
}

int main(int argc, char* argv[])
{
    // Generate basic signal array
    std::vector<std::complex<double>> basic;
    for (unsigned int k= 1; k <= 4; k++) {
        basic.push_back(std::complex<double>(k,0));
    }

    // FFT the thing
    fft(basic);

    for (unsigned int k = 0; k < basic.size(); k++) {
        std::cout << basic.at(k) << std::endl;
    }
}