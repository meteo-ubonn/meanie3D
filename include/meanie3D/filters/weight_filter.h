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

#ifndef M3D_WEIGHTTHRESHOLDFILTER_H
#define M3D_WEIGHTTHRESHOLDFILTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/weights/weight_function.h>
#include <meanie3D/filters/filter.h>

#include <vector>

namespace m3D {

    /** Filters the values in the feature-space based on a set of 
     * thresholds, one for each feature-space dimension. Points that
     * are smaller than the theshold are dropped from the feature-space.
     */
    template<class T>
    class WeightThresholdFilter : public FeatureSpaceFilter<T>
    {
    private:

        WeightFunction <T> *m_weight_function;
        T m_lower_threshold;
        T m_upper_threshold;

    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** Constructor
         * @param weight function
         * @param lower boundary
         * @param upper boundary
         * @param show progress indicator while filtering (default no)
         * @throws logic_error if |thresholds| = 0
         */
        WeightThresholdFilter(WeightFunction <T> *weight_function,
                              T lower_threshold = numeric_limits<T>::min(),
                              T upper_threshold = numeric_limits<T>::max(),
                              bool show_progress = false);

        virtual ~WeightThresholdFilter();

#pragma mark -
#pragma mark Abstract filter method

        virtual void apply(FeatureSpace <T> *fs);
    };
}

#endif
