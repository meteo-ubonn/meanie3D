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

#ifndef M3D_OPERATION_ITERATE_H
#define M3D_OPERATION_ITERATE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/featurespace.h>

namespace m3D { 

    template <typename T>
    class IterationOperation : public Operation<T>
    {

    public:

        IterationOperation( FeatureSpace<T> *fs, PointIndex<T> *index ) : Operation<T>( fs, index ) {};

        virtual ~IterationOperation() {}

        /** Performs a mean-shift iteration from the given starting
         * point until one of the two termination criteria is met.
         * 
         * @param origin of the iteration (a point in feature-space)
         * @param params 
         * @param kernel
         * @param weight weight function to use (or NULL)
         * @param termcrit_epsilon smallest iteration step
         * @param termcrit_iterations maximum number of iterations
         */
        void
        iterate(Point<T> *origin, 
                const SearchParameters *params, 
                const Kernel<T> *kernel, 
                const WeightFunction<T> *weight,
                const T termcrit_epsilon,
                const size_t termcrit_iterations);

        /** Identical to 'iterate' with the difference that the iteration
         * trajectory is recorded and returned.
         * 
         * @param origin of the iteration (a point in feature-space)
         * @param params 
         * @param kernel
         * @param weight weight function to use (or NULL)
         * @param termcrit_epsilon smallest iteration step
         * @param termcrit_iterations maximum number of iterations
         * 
         * @return The trajectory of the iteration
         */
        typename FeatureSpace<T>::Trajectory * 
        get_trajectory(Point<T> *origin, 
                       const SearchParameters *params, 
                       const Kernel<T> *kernel, 
                       const WeightFunction<T> *weight,
                       const T termcrit_epsilon,
                       const size_t termcrit_iterations);


        
    };   
}

#endif
