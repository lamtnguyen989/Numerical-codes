#include <iostream>
#include <complex>
#include <vector>

using Complex = std::complex<float>;

/* Basic recursive FFT (only meant for length of power of 2) */
std::vector<Complex> fft_recurse(std::vector<Complex> X)
{
    // Length and helper vectors
    const unsigned int N = X.size();
    std::vector<Complex> evens(N/2);
    std::vector<Complex> odds(N/2);

    if (N == 1 )
        return X;     // Recursion base case

    // Copy the evens and odds
    for (unsigned int k = 0; k < N/2; k++) {
        evens.at(k) = X.at(2*k);
        odds.at(k) = X.at(2*k + 1);
    }

    // Recurse the DFT
    evens = fft_recurse(evens);
    odds = fft_recurse(odds);

    // Computation result
    for (unsigned int k = 0; k < N/2; k++) {
        Complex root_of_unity = static_cast<Complex>(std::polar(1.0, -2*M_PI*k / N));
        X.at(k) = evens.at(k) + root_of_unity * odds.at(k);
        X.at(k + N/2) = evens.at(k) - root_of_unity * odds.at(k);
    }
    return X;
}

/* Iterative version of above code (power of 2 length) */
std::vector<Complex> fft_iterative_pow_of_2(std::vector<Complex> X)
{
    // Length
    const unsigned int N = X.size();

    // Bit-reversal permutation
    unsigned int index_bit_reversed = 0;
    for (unsigned int index = 0; index < N; index++) {
        // Swapping and make sure we don't double swap
        if (index_bit_reversed > index)
            std::swap(X.at(index), X.at(index_bit_reversed));

        // Calculate the bit reversal of next index
        unsigned int right_shift = N >> 1;
        while ((right_shift >= 1) && (index_bit_reversed >= right_shift)) {
            index_bit_reversed -= right_shift;
            right_shift >> 1;
        }
        index_bit_reversed += right_shift;
    }

    // FFT
    for (unsigned int stage = 2; stage <= N; stage *= 2) {

        // The stage's root of unity
        Complex root_of_unity = static_cast<Complex> (std::polar(1.0, -2*M_PI*stage / N));

        for (unsigned int group = 0; group < N; group += stage) {
            for (unsigned int k = 0; k < stage/2 ; k++) {

                // Calculate the butterfly parts
                Complex p1 = X.at(group + k);
                Complex p2 = root_of_unity * X.at(group + k + stage/2);

                X.at(group + k) =  p1 + p2;     // Lower half
                X.at(group + k + stage/2) = p1 - p2;   // Higher half
            }
        }
    }
    

    return X;
}

/* Wrapper for dealing with arbirary length signals */
std::vector<Complex> fft(std::vector<Complex> X)
{
    // Signal length
    const unsigned int N = X.size();

    // No need to fuss with extra logic if size is a power of 2
    if ((N & (N-1)) == 0)
        return fft_iterative_pow_of_2(X);

    // TODO: Deal with arbitrary length
    return X;
}

/* DFT function to test FFT */
std::vector<Complex> dft(std::vector<Complex> X)
{
    // TODO
}

int main(int argc, char* argv[])
{
    // Generate basic signal array
    std::vector<Complex> basic;
    for (unsigned int k= 1; k <= 4; k++) {
        basic.push_back(Complex(k,0));
    }

    // FFT the thing
    std::vector<Complex> result = fft(basic);

    // Basic (eyes) Testing
    std::cout << "Original signal: " << std::endl;
    for (unsigned int k = 0; k < basic.size(); k++) {
        std::cout << basic.at(k) << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "FFT: " << std::endl;
    for (unsigned int k = 0; k < result.size(); k++) {
        std::cout << result.at(k) << " ";
    }
    std::cout << std::endl;
}
