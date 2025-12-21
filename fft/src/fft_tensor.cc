#include <iostream>
#include <complex>
#include <vector>

using Complex = std::complex<float>;

/* Arbitrary order tensor data structure */
template <typename T>
class Tensor 
{
    public:

        /* Constructors */
        Tensor() {}

        Tensor(std::vector<T> data_values, std::vector<unsigned int> shape_dimensions) 
            : data(data_values)
            , shape(shape_dimensions)
        {
            compute_strides();
        }

        /* Metadata */
        unsigned int order() { return shape.size();}
        unsigned int size() { return data.size();}

        /* Setters */
        void set_value_at(T val, std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            data.at(index) = val;
        }
        
        /* Accessors */
        T value_at(std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            return data.at(index);
        }

    private: 

        /* Data fields */
        std::vector<T> data;
        std::vector<unsigned int> strides;
        std::vector<unsigned int> shape;

        /* Compute strides for constructor */
        void compute_strides() 
        {
            strides.resize(order());
            unsigned int stride = 1;
            for (unsigned int k = 0; k < order; k++) {
                strides.at(k) = stride;
                stride *= shape.at(k);
            }
        }

        /* Process the indices */
        unsigned int process_indices(std::vector<unsigned int>& indices)
        {
            // Checking if the sizes is a-okay
            const unsigned int order = this.order();
            if (indices.size() != order)
                throw std::invalid_argument("Dimension mismatch: Invalid number of indicies provided!");


            // TODO: Check for indices in range 

            // Computing the tensor contiguous data index
            unsigned int index = 0;
            for (unsigned int k = 0; k < order; k++) {
                index += indices.at(k) * strides.at(k);
            }

            return index;
        }
};