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

#ifndef M3D_SCALESPACEKERNEL_IMPL_H
#define M3D_SCALESPACEKERNEL_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include "scalespace_kernel.h"

namespace m3D {

    template<typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const T &t)
            : m_t(t), m_gauging_factor(1.0 / sqrt(2.0 * M_PI * t)), m_values(vector<T>(0)) {
    }

    template<typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const T &t, const vector <T> &distances)
            : m_t(t), m_gauging_factor(1.0 / sqrt(2.0 * M_PI * t)), m_values(vector<T>(distances.size(), 0.0)) {
        for (size_t i = 0; i < distances.size(); i++) {
            m_values[i] = this->value(distances[i]);
        }
    }

    template<typename T>
    ScaleSpaceKernel<T>::ScaleSpaceKernel(const ScaleSpaceKernel<T> &o)
            : m_t(o.m_t), m_gauging_factor(o.m_gauging_factor), m_values(o.m_values) {
    }

    template<typename T>
    ScaleSpaceKernel<T>::~ScaleSpaceKernel() {
    }

    template<typename T>
    ScaleSpaceKernel<T>
    ScaleSpaceKernel<T>::operator=(const ScaleSpaceKernel<T> &other) {
        return ScaleSpaceKernel<T>(other);
    }

#pragma mark -
#pragma mark Accessors

    template<typename T>
    const vector <T> &
    ScaleSpaceKernel<T>::values() {
        return m_values;
    }

    template<typename T>
    bool
    ScaleSpaceKernel<T>::isPreSampled() {
        return m_values.size() > 0;
    }


#pragma mark -
#pragma mark Calculation

    template<typename T>
    const T
    ScaleSpaceKernel<T>::value(size_t index) {
        return m_values[index];
    }

    template<typename T>
    const T
    ScaleSpaceKernel<T>::value(T r) {
        return m_gauging_factor * std::exp(-(r * r) / (2 * m_t));
    }

}

#endif
