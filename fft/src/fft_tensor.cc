#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

#include "tensor.h"

using Complex = std::complex<float>;

/*** 
    1D DFT
***/
std::vector<Complex> dft_1d(std::vector<Complex> &X, bool inverse=false)
{
    const unsigned int N = X.size();

    // Result container
    std::vector<Complex> result(N, Complex(0.0, 0.0));

    for (unsigned int k = 0; k < N; k++) {
        Complex sum = Complex(0.0, 0.0);
        for (unsigned int n = 0; n < N; n++) {

            // Changing the primitive root of unity based on whether or not it is an inverse DFT
            Complex root_of_unity;
            if (inverse == false)
                root_of_unity = static_cast<Complex>(std::polar(1.0, -2*M_PI*k*n / N));
            else
                root_of_unity = static_cast<Complex>(std::polar(1.0, 2*M_PI*k*n / N));


            // DFT sum
            sum += X.at(n) * root_of_unity;
        }

        if (inverse) {sum /= N;}
        result.at(k) = sum;
    }

    return result;
}

// Inverse 1D DFT wrapper
std::vector<Complex> idft_1d(std::vector<Complex> X) { return dft_1d(X, true);}

/*** 
    1D FFT
***/

// Iterative FFT for power of 2 length signals
std::vector<Complex> fft_1d_pow_of_2(std::vector<Complex> X, bool inverse=false)
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
            right_shift >>= 1;
        }
        index_bit_reversed += right_shift;
    }

    // FFT
    for (unsigned int stage = 2; stage <= N; stage *= 2) {

        // The stage's root of unity
        Complex root_of_unity = (inverse == false) ? static_cast<Complex>(std::polar(1.0, -2*M_PI / stage)) 
                                                   : static_cast<Complex>(std::polar(1.0, 2*M_PI / stage));

        for (unsigned int group = 0; group < N; group += stage) {
            Complex twiddle = Complex(1.0, 0.0);
            for (unsigned int k = 0; k < stage/2 ; k++) {

                // Calculate the butterfly parts
                Complex p1 = X.at(group + k);
                Complex p2 = twiddle * X.at(group + k + stage/2);

                X.at(group + k) =  p1 + p2;            // Lower half
                X.at(group + k + stage/2) = p1 - p2;   // Higher half

                // Update twiddle 
                twiddle *= root_of_unity;
            }
        }
    }

    // Inverse transform modification
    if (inverse) {
        for (unsigned int k = 0; k < N; k++) {
            X.at(k) /= N;
        }
    }
    
    return X;
}

// FFT for arbritrary length signals
std::vector<Complex> fft_1d(std::vector<Complex> X, bool inverse=false)
{
    // Signal length
    const unsigned int N = X.size();

    // No need to fuss with extra logic if size is a power of 2
    if ((N & (N-1)) == 0)
        return fft_1d_pow_of_2(X, inverse);

    // TODO: Deal with arbitrary length
    return dft_1d(X, inverse);
}

// Wrapper for inverse FFT
std::vector<Complex> ifft_1d(std::vector<Complex> X) { return fft_1d(X, true);}


/***
    Arbitrary order DFT using Tensor
***/
// Global indices grabbing
std::vector<unsigned int> grab_tensor_index_of_slice(unsigned int dim, unsigned int slice, unsigned int order, std::vector<unsigned int>& shape)
{
        // Indices container
        std::vector<unsigned int> indx(order, 0);
        
        // Calculate indices for this specific slice
        unsigned int temp = slice;
        for (unsigned int n = 0; n < order; n++) {
                if (n == dim)
                    continue;
                unsigned int stride_val = 1;
                for (unsigned int k = n + 1; k < order; k++) {
                        if (k != dim) 
                            stride_val *= shape.at(k);
                }
                indx.at(n) = (temp / stride_val) % shape.at(n);
        } 
        
        return indx;
}

// Forward transformation
Tensor<Complex> dft(Tensor<Complex> X, bool inverse=false)
{
        // Tensor metadata
        std::vector<unsigned int> shape = X.shape();
        unsigned int order = X.order();
        unsigned int total_size = X.size();
        
        // Iteratively apply 1D DFT along the dimensions
        for (unsigned int dim = 0; dim < order; dim++) {
                // Dimension size and num of slices needed
                unsigned int dim_size = shape.at(dim);
                unsigned int num_slices = total_size / dim_size;
                
                // Result container for this dimension
                Tensor<Complex> result(shape);
                
                // DFT through all of the slices
                for (unsigned int slice = 0; slice < num_slices; slice++) {
                    std::vector<unsigned int> indices = grab_tensor_index_of_slice(dim, slice, order, shape);
                    std::vector<Complex> slice_data = X.extract_1d_slice(dim, indices);
                    std::vector<Complex> transformed = dft_1d(slice_data, inverse);
                    result.set_1d_slice(dim, indices, transformed);
                }
                
                X = result;
        }
        
        return X;
}

// Inverse transformation
Tensor<Complex> idft(Tensor<Complex> X) { return dft(X, true);}


/***
    Arbitrary order FFT using Tensor
***/
// Forward transformation
Tensor<Complex> fft(Tensor<Complex> X, bool inverse=false)
{
        // Tensor metadata
        std::vector<unsigned int> shape = X.shape();
        unsigned int order = X.order();
        unsigned int total_size = X.size();
        
        // Iteratively apply 1D DFT along the dimensions
        for (unsigned int dim = 0; dim < order; dim++) {
                // Dimension size and num of slices needed
                unsigned int dim_size = shape.at(dim);
                unsigned int num_slices = total_size / dim_size;
                
                // Result container for this dimension
                Tensor<Complex> result(shape);
                
                // DFT through all of the slices
                for (unsigned int slice = 0; slice < num_slices; slice++) {
                    std::vector<unsigned int> indices = grab_tensor_index_of_slice(dim, slice, order, shape);
                    std::vector<Complex> slice_data = X.extract_1d_slice(dim, indices);
                    std::vector<Complex> transformed = fft_1d(slice_data, inverse);
                    result.set_1d_slice(dim, indices, transformed);
                }
                
                X = result;
        }
        
        return X;
}


// Inverse transformation
Tensor<Complex> ifft(Tensor<Complex> X) { return fft(X, true);}



int main()
{
    Tensor<Complex> tensor({Complex(1.0, 2.0), Complex(3.0, 4.0), Complex(5.0, 6.0), Complex(7.0, 8.0)}, {2,2});
    std::cout << "Original Tensor:" << std::endl;
    tensor.print_tensor();

    Tensor<Complex> dft_result = dft(tensor);
    Tensor<Complex> fft_result = fft(tensor);

    std::cout << "DFT result: " << std::endl;
    dft_result.print_tensor();

    std::cout << "FFT result: " << std::endl;
    fft_result.print_tensor();

    dft_result = dft(dft_result, true);
    std::cout << "Inverse DFT back result: " << std::endl;
    dft_result.print_tensor();

    fft_result = ifft(fft_result);
    std::cout << "Inverse FFT back result: " << std::endl;
    fft_result.print_tensor();
}