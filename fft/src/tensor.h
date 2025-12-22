#ifndef TENSOR_H
#define TENSOR_H

#include <iostream>
#include <vector>


/* Arbitrary order tensor data structure */
template <typename T>
class Tensor 
{
    public:

        /***
            Constructors 
        ***/
        Tensor() {}

        Tensor(std::vector<T> data_values, std::vector<unsigned int> shape_dimensions) 
            : _data(data_values)
            , _shape(shape_dimensions)
        {
            if (size_from_shape(shape_dimensions) != data_values.size())
                throw std::invalid_argument("Dimension mismatch!");
            
            compute_strides();
        }

        Tensor(std::vector<unsigned int> shape)
            : _shape(shape)
        {
            unsigned int total_size = size_from_shape(shape);
            _data.resize(total_size, T() );
            compute_strides();
        }

        /*** 
            Metadata 
        ***/
        unsigned int order() { return _shape.size();}
        unsigned int size() { return _data.size();}
        std::vector<unsigned int> shape() { return _shape;}

        /*** 
            Setters 
        ***/
        void set_value_at(T val, std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            _data.at(index) = val;
        }

        void set_1d_slice(unsigned int dim, std::vector<unsigned int> indices, std::vector<T>& new_values) 
        {
            unsigned int slice_size = _shape.at(dim);

            for (unsigned int k = 0; k < slice_size; k++) {
                indices.at(dim) = k;
                set_value_at(new_values.at(k), indices);
            }
        }
        
        /*** 
            Accessors 
        ***/
        T value_at(std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            return _data.at(index);
        }

        std::vector<T> extract_1d_slice(unsigned int dim, std::vector<unsigned int> indices)
        {
            unsigned int slice_size = _shape.at(dim);
            std::vector<T> slice(slice_size);

            for (unsigned int k = 0; k < slice_size; k++) {
                indices.at(dim) = k;
                slice.at(k) = value_at(indices);
            }

            return slice;
        }

    private: 

        /***
            Data fields 
        ***/
        std::vector<T> _data;
        std::vector<unsigned int> _strides;
        std::vector<unsigned int> _shape;

        /*** 
            Constructor helpers
        ***/
        void compute_strides() 
        {
            // Configure the strides data container
            unsigned int order = _shape.size();
            _strides.resize(order);

            // Stride with Row-major access
            unsigned int stride_val = 1;
            for (int k = order - 1; k >= 0; k--) {
                _strides.at(k) = stride_val;
                stride_val *= _shape.at(k);
            }
        }

        unsigned int size_from_shape(std::vector<unsigned int>& shape)
        {
            unsigned int total_size = 1;

            for (unsigned int k = 0; k < shape.size(); k++) {
                total_size *= shape.at(k);
            }

            return total_size;
        }

        /*** 
            Process the indices 
        ***/
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