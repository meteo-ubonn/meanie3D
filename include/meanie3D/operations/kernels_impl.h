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

#ifndef M3D_KERNELS_IMPL_H
#define M3D_KERNELS_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils/vector_utils.h>

#include <hdf5.h>
#include <iostream>

#include "kernels.h"

namespace m3D {

    using utils::vectors::vector_norm;

    // Gaussian

    template<class T>
    T GaussianNormalKernel<T>::apply(const T dist) const {
        return 0.5 * exp(-0.5 * dist);
    }

    template<class T>
    T GaussianNormalKernel<T>::apply(const vector <T> &c) const {
        return this->apply(vector_norm(c));
    }

    // Epanechnikov

    template<class T>
    T EpanechnikovKernel<T>::apply(const T dist) const {
        return (dist <= this->m_kernelSize ? this->m_kernelSize - dist : 0);
    }

    template<class T>
    T EpanechnikovKernel<T>::apply(const vector <T> &c) const {
        return this->apply(vector_norm(c));
    }

    // Uniform

    template<class T>
    T UniformKernel<T>::apply(const T dist) const {
        return (fabs(dist) <= this->m_kernelSize) ? 1.0f / this->m_kernelSize : 0;
    }

    template<class T>
    T UniformKernel<T>::apply(const vector <T> &c) const {
        return this->apply(vector_norm(c));
    }
}

#endif
