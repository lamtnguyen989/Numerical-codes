#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <vector>


/* Arbitrary order tensor data structure */
template <typename T>
class Tensor 
{
    public:

        /* Constructors */
        Tensor() {}

        Tensor(std::vector<T> data_values, std::vector<unsigned int> shape_dimensions) 
            : _data(data_values)
            , _shape(shape_dimensions)
        {
            compute_strides();
        }

        /* Metadata */
        unsigned int order() { return _shape.size();}
        unsigned int size() { return _data.size();}
        std::vector<unsigned int> shape() { return _shape;}

        /* Setters */
        void set_value_at(T val, std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            _data.at(index) = val;
        }
        
        /* Accessors */
        T value_at(std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            return _data.at(index);
        }

    private: 

        /* Data fields */
        std::vector<T> _data;
        std::vector<unsigned int> _strides;
        std::vector<unsigned int> _shape;

        /* Compute strides for constructor */
        void compute_strides() 
        {
            // Configure the strides data container
            unsigned int order = _shape.size();
            _strides.resize(order);

            // Stride with Row-major access
            _strides.back() = 1;
            if ((int) order - 2 >= 0) {
                for (unsigned int k = order - 2; k >= 0; k--) {
                    _strides.at(k) = _strides.at(k+1) * _shape.at(k+1);
                }
            }
        }

        /* Process the indices */
        unsigned int process_indices(std::vector<unsigned int>& indices)
        {
            // Checking if the sizes is a-okay
            const unsigned int order = _shape.size();
            if (indices.size() != order)
                throw std::invalid_argument("Dimension mismatch!");

            // Computing the tensor contiguous data index
            unsigned int index = 0;
            for (unsigned int k = 0; k < order; k++) {

                // Checking every index is within bounds
                if (indices.at(k) >= _shape.at(k))
                    throw std::out_of_range("Index out of bounds!");

                // Update the global flat index
                index += indices.at(k) * _strides.at(k);
            }

            return index;
        }
};

#endif