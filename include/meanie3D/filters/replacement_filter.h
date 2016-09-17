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

#ifndef M3D_REPLACEMENT_FILTER_H
#define M3D_REPLACEMENT_FILTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include "filter.h"

namespace m3D {

    /**
     * This filter replaces the value at each point in featurspace with a value
     * calculated from it's neighbours. It takes generally four parameters. One
     * is the mode, which can take those three values:
     * <ul>
     * <li>ReplaceWithLowest - replaces the value with the average of the x% lowest neighbours</li>
     * <li>ReplaceWithHighest - replaces the value with the average of the x% highest neighbours</li>
     * <li>ReplaceWithMedian - replaces the value with the median of all neighbours</li>
     * </ul>
     *
     * Another parameter is the index of variable in the featurespace's value range that
     * is supposed to be processed.
     *
     * Next the bandwidth (the radius of the filter). The last is a percentage which is
     * used in either ReplaceWithLowest or ReplaceWithHighest (see above).
     *
     * Note that the arithmetic average can be calculated by using either lowest
     * or highest and 100%.
     */
    template<class T>
    class ReplacementFilter : public FeatureSpaceFilter<T>
    {
    public:

        typedef enum
        {
            ReplaceWithLowest,
            ReplaceWithHighest,
            ReplaceWithMedian
        } ReplacementMode;

    private:

        ReplacementMode m_replacement_mode;
        size_t m_variable_index;
        float m_percentage;
        std::vector<T> m_bandwidth;

    public:

        ReplacementFilter(const ReplacementMode mode,
                          const size_t variable_index,
                          const std::vector<T> &bandwidth,
                          const float percentage = 0.25,
                          bool show_progress = false);

#pragma mark -
#pragma mark Constructor/Destructor

        virtual ~ReplacementFilter() {}

#pragma mark -
#pragma mark Abstract filter method

        virtual void apply(FeatureSpace<T> *fs);
    };
}

#endif
