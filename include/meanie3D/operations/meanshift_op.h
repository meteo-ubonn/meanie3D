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

#ifndef M3D_OPERATION_MEANSHIFTOPERATION_H
#define M3D_OPERATION_MEANSHIFTOPERATION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/operations.h>
#include <meanie3D/featurespace.h>
#include <meanie3D/index.h>

#include <vector>

namespace m3D {

    template <typename T>
    class MeanshiftOperation : public Operation<T>
    {
    public:

        //static const int NO_WEIGHT;

        MeanshiftOperation(FeatureSpace<T> *fs,
                PointIndex<T> *index) : Operation<T>(fs, index)
        {
        }

        virtual ~MeanshiftOperation()
        {
        }


        /** Call this up-front to preempt lazy index construction.
         * @param search parameters
         */
        void prime_index(const SearchParameters *params);

        /** Meanshift calculation at point x. Sample is done by range search. 
         * {@see SampleOperation::sample_range}
         * 
         * @param feature space coordinate x (origin)
         * @param search parameters
         * @param kernel
         * @param weight function (defaults to NULL)
         * @param flag, indicating if the returned vector should be rounded
         *        to the coordinate system's resolution
         * @return mean shift vector
         */
        vector<T>
        meanshift(const vector<T> &x,
                const SearchParameters *params,
                const Kernel<T> *kernel = new GaussianNormalKernel<T>(),
                const WeightFunction<T> *w = NULL,
                const bool normalize_shift = true);
    };
}

#endif
