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
            : data_(data_values)
            , shape_(shape_dimensions)
        {
            compute_strides();
        }

        /* Metadata */
        unsigned int order() { return shape_.size();}
        unsigned int size() { return data_.size();}
        std::vector<unsigned int> shape() { return shape_;}

        /* Setters */
        void set_value_at(T val, std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            data_.at(index) = val;
        }
        
        /* Accessors */
        T value_at(std::vector<unsigned int>& indices)
        {
            unsigned int index = process_indices(indices);
            return data_.at(index);
        }

    private: 

        /* Data fields */
        std::vector<T> data_;
        std::vector<unsigned int> strides_;
        std::vector<unsigned int> shape_;

        /* Compute strides for constructor */
        void compute_strides() 
        {
            // Configure the strides data container
            unsigned int order = shape_.size();
            strides_.resize(order);

            // Stride with Row-major access
            strides_.back() = 1;
            for (unsigned int k = order - 2; k >= 0; k--) {
                strides_.at(k) = strides_.at(k+1) * shape_.at(k+1);
            }
        }

        /* Process the indices */
        unsigned int process_indices(std::vector<unsigned int>& indices)
        {
            // Checking if the sizes is a-okay
            const unsigned int order = this.order();
            if (indices.size() != order)
                throw std::invalid_argument("Dimension mismatch: Indicies size must aligned with the order of tensors!");

            // Computing the tensor contiguous data index
            unsigned int index = 0;
            for (unsigned int k = 0; k < order; k++) {

                // Checking every index is within bounds
                if (indices.at(k) >= shape_.at(k))
                    throw std::out_of_range("Index out of bounds!");

                // Update the global flat index
                index += indices.at(k) * strides_.at(k);
            }

            return index;
        }
};