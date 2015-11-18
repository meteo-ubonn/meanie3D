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

#ifndef M3D_OPERATION_MEANSHIFTOPERATION_IMPL_H
#define M3D_OPERATION_MEANSHIFTOPERATION_IMPL_H

#include <stdexcept>
#include <cmath>
#include <netcdf>

#include "meanshift_op.h"

namespace m3D {

    template <typename T>
    void
    MeanshiftOperation<T>::prime_index(const SearchParameters *params)
    {
        vector<T> x(this->feature_space->rank());
        typename Point<T>::list *sample = this->point_index->search(x, params, NULL);
        delete sample;
    }

    template <typename T>
    vector<T>
    MeanshiftOperation<T>::meanshift(const vector<T> &x,
            const SearchParameters *params,
            const Kernel<T> *kernel,
            const WeightFunction<T> *w,
            const bool normalize_shift)
    {
        using namespace utils::vectors;

        // Distances are only used in the KNN case.
        vector<T> distances;
        typename Point<T>::list *sample = NULL;
        if (params->search_type() == SearchTypeRange) {
            sample = this->point_index->search(x, params, NULL);
        } else {
            sample = this->point_index->search(x, params, &distances);
        }
        vector<T> shift(this->feature_space->dimension, 0.0);

        // If the sample is empty, no shift can be calculated.
        // Returns a shift of 0
        if (sample->size() == 0) {
            return shift;
        }

        vector<T> numerator(this->feature_space->dimension);
        T denominator = 0.0;
        size_t size = sample->size();

        // Range Search
        if (params->search_type() == SearchTypeRange) {
            vector<T> h = ((RangeSearchParams<T> *) params)->bandwidth;
            if (w == NULL && kernel != NULL) {
                // No weight function
                for (size_t index = 0; index < size; index++) {
                    T weight = kernel->apply(mahalabonis_distance_sqr(x, sample->at(index)->values, h));
                    denominator += weight;
                    numerator += weight * sample->at(index)->values;
                }
            } else if (w != NULL && kernel == NULL) {
                // no kernel
                for (size_t index = 0; index < size; index++) {
                    T weight = w->operator()(sample->at(index));
                    denominator += weight;
                    numerator += weight * sample->at(index)->values;
                }
            } else if (kernel == NULL && w == NULL) {
                // no weight function AND no kernel
                for (size_t index = 0; index < size; index++) {
                    denominator += 1.0;
                    numerator += sample->at(index)->values;
                }
            } else {
                // Both weight function and kernel
                vector<T> mult(x.size(), 0.0);
                for (size_t index = 0; index < size; index++) {
                    T var_weight = w->operator()(sample->at(index));
                    T kernel_weight = kernel->apply(mahalabonis_distance_sqr(x, sample->at(index)->values, h));
                    T weight = kernel_weight * var_weight;
                    denominator += weight;
                    for (int i = 0; i < x.size(); i++) {
                        numerator[i] += (weight * sample->at(index)->values[i]);
                    }
                }
            }
        }            // KNN
        else {
            throw "Not Implemented";
        }

        vector<T> dx = numerator / denominator;
        shift = dx - x;
        delete sample;
        if (normalize_shift) {
            shift = this->feature_space->coordinate_system->round_to_grid(shift);
        }

        return shift;
    }
}

#endif
