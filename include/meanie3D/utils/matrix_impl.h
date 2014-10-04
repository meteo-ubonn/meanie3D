#ifndef _M3D_Matrix_Impl_H_
#define _M3D_Matrix_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include "matrix.h"

namespace m3D { namespace utils {

    template <typename T>
    typename Matrix<T>::matrix_t
    Matrix<T>::create_matrix(size_t width, size_t height)
    {
        matrix_t matrix;

        matrix.resize(width);

        for (int i=0; i<width; ++i)
        {
            matrix[i].resize(height);
        }

        return matrix;
    }

    template <typename T>
    typename Matrix<T>::flag_matrix_t
    Matrix<T>::create_flag_matrix(size_t width, size_t height, int defaultValue)
    {
        flag_matrix_t matrix;

        matrix.resize(width);

        for (size_t i=0; i<width; ++i)
        {
            matrix[i].resize(height,defaultValue);
        }

        return matrix;
    }
}}

#endif
