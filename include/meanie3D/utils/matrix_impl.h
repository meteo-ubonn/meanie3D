/* The MIT License (MIT)
 * 
 * (c) Jürgen Simon 2014 (juergen.simon@uni-bonn.de)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _M3D_Matrix_Impl_H_
#define _M3D_Matrix_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include "matrix.h"

namespace m3D { namespace utils {

    template <typename T>
    typename SimpleMatrix<T>::matrix_t
    SimpleMatrix<T>::create_matrix(size_t width, size_t height)
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
    typename SimpleMatrix<T>::flag_matrix_t
    SimpleMatrix<T>::create_flag_matrix(size_t width, size_t height, int defaultValue)
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
