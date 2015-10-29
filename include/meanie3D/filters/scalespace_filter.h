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

#ifndef M3D_SCALESPACEFILTER_H
#define M3D_SCALESPACEFILTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/filters/filter.h>
#include <meanie3D/filters/scalespace_kernel.h>

#include <boost/progress.hpp>

namespace m3D {

    /** Smoothes the data with a scale-space filter. This filter does NOT create
     * new points. Only the existing points are smoothed out.
     */
    template <typename T>
    class ScaleSpaceFilter : public FeatureSpaceFilter<T>
    {
    private:

        T m_scale; /// Scale

        T m_decay; /// Decay

        vector<ScaleSpaceKernel<T> > m_kernels; /// Scale-Space kernel

        map<size_t, T> m_unfiltered_min; /// unfiltered minimum
        map<size_t, T> m_unfiltered_max; /// unfiltered maximum

        map<size_t, T> m_min; /// minimum tracker
        map<size_t, T> m_max; /// maximum tracker

        boost::progress_display *m_progress_bar; /// Progress meter

        size_t m_modified_points;

        size_t m_created_points;

        vector<netCDF::NcVar> m_excluded_vars;

        void
        applyWithArrayIndexRecursive(FeatureSpace<T> *fs,
                ArrayIndex<T> *originalIndex,
                ArrayIndex<T> *filteredPoints,
                vector<size_t> &dimensionIndexes,
                size_t dimensionIndex,
                typename CoordinateSystem<T>::GridPoint& gridpoint);

        void
        applyWithArrayIndexForDimension(FeatureSpace<T> *fs,
                ArrayIndex<T> *originalIndex,
                ArrayIndex<T> *filteredPoints,
                size_t fixedDimensionIndex);

        void
        applyWithArrayIndex(FeatureSpace<T> *fs);

#pragma mark -
#pragma mark Paralellized versions

        void
        apply_parallellized_on_dimension(FeatureSpace<T> *fs,
                ArrayIndex<T> *originalIndex,
                ArrayIndex<T> *filteredPoints,
                size_t fixedDimensionIndex);

        void
        apply_parallellized(FeatureSpace<T> *fs);

    public:

#pragma mark -
#pragma mark Constructor/Destructor

        /** Constructor
         * @param scale parameter
         * @param resolution vector
         * @param decay in percent of value at 0
         * @param show progress indicator while filtering (default no)
         * @throws logic_error if scale < 0 or decay < 0 or > 1.
         */
        ScaleSpaceFilter(T scale,
                const vector<T> &resolution,
                vector<netCDF::NcVar> &exclude_from_scale_space_filtering,
                T decay = 0.01,
                bool show_progress = false);

        /** Destructor
         */
        virtual ~ScaleSpaceFilter();

#pragma mark -
#pragma mark Abstract filter method

        /** Applies the scale-space filter. Modifies the given
         * instance of F.
         *
         * @param feature space
         */
        virtual void apply(FeatureSpace<T> *fs);

#pragma mark -
#pragma mark After the processing 

        /** After the filtering, one can retrieve the value ranges
         * of the filtered scan from these variables. 
         * @return filtered minimum
         */
        const map<size_t, T> &
        get_filtered_min();

        /** After the filtering, one can retrieve the value ranges
         * of the filtered scan from these variables.
         * @return filtered maximum
         */
        const map<size_t, T> &
        get_filtered_max();

        /** Retrieve the scaling ratios for the individual variables
         * after the filter has run. Each ratio is calculated as
         * 
         *      ratio = filtered range / unfiltered range 
         * 
         * This describes how the range is changed under transformation.
         */
        map<size_t, T>
        getRangeFactors();

    };
}

#endif
