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

    template <class T>
    T GaussianNormalKernel<T>::apply(const T dist) const {
        return 0.5 * exp(-0.5 * dist);
    }

    template <class T>
    T GaussianNormalKernel<T>::apply(const vector<T> &c) const {
        return this->apply(vector_norm(c));
    }

    // Epanechnikov

    template <class T>
    T EpanechnikovKernel<T>::apply(const T dist) const {
        return ( dist <= this->m_kernelSize ? this->m_kernelSize - dist : 0);
    }

    template <class T>
    T EpanechnikovKernel<T>::apply(const vector<T> &c) const {
        return this->apply(vector_norm(c));
    }

    // Uniform

    template <class T>
    T UniformKernel<T>::apply(const T dist) const {
        return ( fabs(dist) <= this->m_kernelSize) ? 1.0f / this->m_kernelSize : 0;
    }

    template <class T>
    T UniformKernel<T>::apply(const vector<T> &c) const {
        return this->apply(vector_norm(c));
    }
}

#endif
