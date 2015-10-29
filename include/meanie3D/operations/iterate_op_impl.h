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

#ifndef M3D_OPERATION_ITERATE_IMPL_H
#define M3D_OPERATION_ITERATE_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include "iterate_op.h"

namespace m3D {

#pragma mark -
#pragma mark Meanshift Iteration    

    template <typename T>
    typename FeatureSpace<T>::Trajectory *
    IterationOperation<T>::get_trajectory(Point<T> *origin,
            const SearchParameters *params,
            const Kernel<T> *kernel,
            const WeightFunction<T> *weight,
            const T termcrit_epsilon,
            const size_t termcrit_iterations)
    {
        using namespace m3D::vectors;

        // Allocate a fresh trajectory

        typename FeatureSpace<T>::Trajectory *trajectory = new typename FeatureSpace<T>::Trajectory();

        // Add the origin as first point

        vector<T> x = origin->values;

        trajectory->push_back(x);

        // Go

        T dx = std::numeric_limits<T>::max();

        vector<T> d(x.size());

        size_t iter = 0;

        while (iter < termcrit_iterations && dx >= termcrit_epsilon) {
#if DEBUG_ITERATION
            std::cout << "Iteration " << iter << " from " << x;
#endif
            // get the mean-shift

            MeanshiftOperation<T> op(this->feature_space, this->point_index);

            vector<T> shift = op.meanshift(x, params, kernel, weight);

            if (iter == 0) {
                origin->shift = shift;
            }

#if DEBUG_ITERATION
            std::cout << " with mean-shift " << shift;
#endif            

            // calculate termination criteria
            dx = (T) vector_norm(shift);

            // calculate iteration end point (re-use shift variable)
            // Use a slightly optimized loop form 

            size_t index = 0;
            typename vector<T>::iterator it;
            for (it = shift.begin(); it != shift.end(); it++) {
                *it += x[index++];
            }

#if DEBUG_ITERATION
            std::cout << " to " << shift;
#endif            
            // add to trajectory

            trajectory->push_back(shift);

#if DEBUG_ITERATION
            std::cout << " (dx = " << dx << ")" << std::endl;
#endif            
            iter++;

            x = shift;
        }

        // substract one, since termination criterion is checked
        // at beginning of the loop
        iter--;

        return trajectory;
    }

    template <typename T>
    void
    IterationOperation<T>::iterate(Point<T> *origin,
            const SearchParameters *params,
            const Kernel<T> *kernel,
            const WeightFunction<T> *weight,
            const T termcrit_epsilon,
            const size_t termcrit_iterations)
    {
        using namespace m3D::vectors;

        vector<T> x = origin->values;

        T dx = std::numeric_limits<T>::max();

        size_t iter = 0;

        while (iter < termcrit_iterations && dx >= termcrit_epsilon) {
#if DEBUG_ITERATION
            std::cout << "Iteration " << iter << " from " << x;
#endif
            // get the mean-shift

            MeanshiftOperation<T> op(this->feature_space, this->point_index);

            vector<T> shift = op.meanshift(x, params, kernel, weight);

            if (iter == 0) {
                origin->shift = shift;
            }

#if DEBUG_ITERATION
            std::cout << " with mean-shift " << shift;
#endif            
            // calculate termination criteria

            dx = (T) vector_norm(shift);

            // calculate iteration end point (re-use shift variable)
            // Use a slightly optimized loop form 

            size_t index = 0;
            typename vector<T>::iterator it;
            for (it = shift.begin(); it != shift.end(); it++) {
                *it += x[index++];
            }

#if DEBUG_ITERATION
            std::cout << " to " << shift;
            std::cout << " (dx = " << dx << ")" << std::endl;
#endif            
            iter++;

            x = shift;
        }

        // substract one, since termination criterion is checked
        // at beginning of the loop
        iter--;
    }
}

#endif
