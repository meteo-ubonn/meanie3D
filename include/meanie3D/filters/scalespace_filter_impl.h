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

#ifndef M3D_SCALESPACEFILTER_IMPL_H
#define M3D_SCALESPACEFILTER_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/utils.h>

#include <exception>
#include <stdexcept>
#include <cmath>
#include <map>
#include <string>

#include <boost/progress.hpp>

#include "scalespace_filter.h"

namespace m3D {

    using namespace std;

#pragma mark -
#pragma mark Utility methods

    template <typename T>
    T ScaleSpaceFilter<T>::scale_to_filter_width(T t, T decay) {
        T filter_width = sqrt(ceil(-2.0 * t * log(decay))) / 2.0;
        return filter_width;
    }

    template <typename T>
    T ScaleSpaceFilter<T>::filter_width_to_scale(T width, T decay) {
        T scale = -2.0 * width * width / log(decay);
        return scale;
    }

#pragma mark -
#pragma mark Constructor/Destructor

    template <typename T>
    ScaleSpaceFilter<T>::ScaleSpaceFilter(T scale,
            const vector<T> &resolution,
            const vector<NcVar> &excluded_vars,
            T decay,
            bool show_progress)
    : FeatureSpaceFilter<T>(show_progress)
    , m_scale(scale)
    , m_decay(decay)
    , m_progress_bar(NULL)
    , m_excluded_vars(excluded_vars)
    {
        if (scale < 0) {
            throw logic_error("scale can not be less than zero");
        }

        if (decay < 0 || decay >= 1) {
            throw logic_error("decay must be > 0 and < 1");
        }

        T filter_width = sqrt(ceil(-2.0 * scale * log(decay))) / 2;

        m_kernels.clear();

        for (size_t i = 0; i < resolution.size(); i++) {
            // calculate the distances vector

            size_t mask_size = filter_width / resolution[i];

            vector<T> distances(mask_size, 0.0);

            for (size_t j = 0; j < mask_size; j++) {
                distances[j] = ((T) j) * resolution[i];
            }

            // create the kernel

            ScaleSpaceKernel<T> kernel(scale, distances);

            m_kernels.push_back(kernel);
        }
    }

    template <typename T>
    ScaleSpaceFilter<T>::~ScaleSpaceFilter()
    {
    }

#pragma mark -
#pragma mark Abstract filter method

    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndexRecursive(FeatureSpace<T> *fs,
            ArrayIndex<T> *originalIndex,
            ArrayIndex<T> *filteredPoints,
            vector<size_t> &dimensionIndexes,
            size_t dimensionIndex,
            typename CoordinateSystem<T>::GridPoint& gridpoint)
    {
        // loop over the dimension with the given index. For each
        // coordinate, loop over the other two. At Each point, apply
        // the one-dimensional kernel

        using namespace std;

        const CoordinateSystem<T> *cs = fs->coordinate_system;

        size_t realDimIndex = dimensionIndexes[dimensionIndex];

        NcDim dim = cs->dimensions()[realDimIndex];

        // iterate over dimensions

        for (int index = 0; index < dim.getSize(); index++) {
            gridpoint[realDimIndex] = index;

            if (dimensionIndex < (gridpoint.size() - 1)) {
                applyWithArrayIndexRecursive(fs, originalIndex, filteredPoints, dimensionIndexes, dimensionIndex + 1, gridpoint);
            } else {
                // we reached the fixed dimension

                if (this->show_progress()) {
                    m_progress_bar->operator++();
                }

                // exclude points that were off limits
                // in any of the original data sets

                //if (fs->off_limits()->get(gridpoint))
                //    continue;

                ScaleSpaceKernel<T> g = this->m_kernels[realDimIndex];

                // Find the boundaries. Take care not to step
                // outside the bounds of the array

                int width = g.values().size() - 1;
                int gpIndex = (int) gridpoint[realDimIndex];
                int minIndex = (gpIndex - width >= 0) ? (gpIndex - width) : 0;
                int maxIndex = ((gpIndex + width) < (dim.getSize() - 1)) ? (gpIndex + width) : (dim.getSize() - 1);

                // Convolute in 1D around the given point with
                // the mask size determined by the kernel
                // Run the convolution for each feature variable

                typename CoordinateSystem<T>::GridPoint gridIter = gridpoint;

                vector<T> sum(fs->value_rank(), 0.0);

                size_t sumCount = 0;

                for (int i = minIndex; i < maxIndex; i++) {
                    // set gridpoint to current position

                    gridIter[realDimIndex] = i;

                    // Again, make sure no points originally marked as
                    // off limits are used

                    //                    if (fs->off_limits()->get(gridIter))
                    //                        continue;

                    // get the point at the iterated position

                    Point<T> *pIter = originalIndex->get(gridIter);

                    if (pIter == NULL)
                        continue;

                    // apply the pre-sampled gaussian and sum up

                    size_t d = (i <= index) ? (index - i) : (i - index);

                    for (int varIndex = 0; varIndex < fs->value_rank(); varIndex++) {
                        T value = pIter->values[cs->rank() + varIndex];
                        sum[varIndex] += g.value(d) * value;
                    }

                    sumCount++;
                }

                // No muss, no fuss

                if (sumCount == 0)
                    continue;

                // Fuss! Fetch the point from the array index

                Point<T> *p = filteredPoints->get(gridpoint);

                // If no point existed, decide if we need to create one

                if (p == NULL) {
                    // Create a new point with default values
                    // and insert into array index

                    typename CoordinateSystem<T>::Coordinate coordinate = cs->newCoordinate();
                    cs->lookup(gridpoint, coordinate);
                    vector<T> values = coordinate;
                    values.resize(fs->rank(), 0.0);

                    p = PointFactory<T>::get_instance()->create(gridpoint, coordinate, values);

                    // Did this exist in the original index?

                    Point<T> *op = originalIndex->get(gridpoint);
                    p->isOriginalPoint = ((op == NULL) ? false : op->isOriginalPoint);

                    // Since we just created this point, there
                    // is no need to copy it again when adding
                    // it to the array index

                    filteredPoints->set(gridpoint, p, false);
                }

                // If we have a point after all that, update it with the
                // filtered value

                if (p != NULL) {
                    // copy values and track limits

                    for (int varIndex = 0; varIndex < fs->value_rank(); varIndex++) {
                        p->values[fs->spatial_rank() + varIndex] = sum[varIndex];

                        if (sum[varIndex] < m_min[varIndex])
                            m_min[varIndex] = sum[varIndex];

                        if (sum[varIndex] > m_max[varIndex])
                            m_max[varIndex] = sum[varIndex];
                    }
                }
            }
        }
    }

    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndexForDimension(FeatureSpace<T> *fs,
            ArrayIndex<T> *originalIndex,
            ArrayIndex<T> *filteredPoints,
            size_t fixedDimensionIndex)
    {
        // iterate over the other two

        // Create a vector, that enumerates the dimensions such, that the
        // fixed dimension comes last

        vector<size_t> dimensionIndexes;

        for (size_t j = 0; j < fs->coordinate_system->rank(); j++) {
            if (j == fixedDimensionIndex) continue;
            dimensionIndexes.push_back(j);
        }

        dimensionIndexes.push_back(fixedDimensionIndex);

        typename CoordinateSystem<T>::GridPoint gridpoint = fs->coordinate_system->newGridPoint();

        // Now recurse into the structure, bearing in mind that
        // the dimensions have been re-ordered

        applyWithArrayIndexRecursive(fs, originalIndex, filteredPoints, dimensionIndexes, 0, gridpoint);
    }

    template <typename T>
    void
    ScaleSpaceFilter<T>::applyWithArrayIndex(FeatureSpace<T> *fs)
    {
        using namespace std;

        const CoordinateSystem<T> *cs = fs->coordinate_system;

        // index the original

        if (this->show_progress()) {
            cout << endl << "Constructing array indexes ...";

            start_timer();
        }

        ArrayIndex<T> *originalIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), fs->points, true);

        ArrayIndex<T> *filteredIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), false);

        if (this->show_progress()) {
            cout << "done. (" << stop_timer() << "s)" << endl;
        }

        if (this->show_progress()) {
            cout << endl << "Applying scale filter t=" << m_scale << " decay=" << m_decay << " ... " << endl;

            long numPoints = 1;

            for (size_t i = 0; i < fs->coordinate_system->rank(); i++) {
                numPoints *= fs->coordinate_system->dimensions()[i].getSize();
            }

            m_progress_bar = new boost::progress_display(fs->spatial_rank() * numPoints);

            start_timer();
        }

        // initialize min/max and re-set counts

        for (size_t varIndex = 0; varIndex < fs->value_rank(); varIndex++) {
            m_min[varIndex] = std::numeric_limits<T>::max();
            m_max[varIndex] = std::numeric_limits<T>::min();
        }

        // Apply dimension by dimension (exploiting separability)

        for (size_t dimIndex = 0; dimIndex < fs->spatial_rank(); dimIndex++) {
            applyWithArrayIndexForDimension(fs, originalIndex, filteredIndex, dimIndex);

            delete originalIndex;

            if (dimIndex < (fs->spatial_rank() - 1)) {
                originalIndex = filteredIndex;

                filteredIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), false);
            }
        }

        // replace the points in the original with the filtered
        // array index results
        filteredIndex->replace_points(fs->points);

        size_t originalPoints = 0;
        for (size_t i = 0; i < fs->points.size(); i++) {
            if (fs->points[i]->isOriginalPoint) originalPoints++;
        }

        if (this->show_progress()) {
            cout << "done. (" << stop_timer() << "s)" << endl;
            cout << "Filtered featurespace contains " << fs->size() << " points (" << originalPoints << " original points, "
                    << "(" << (fs->size() - originalPoints) << " new points))" << endl;
            delete m_progress_bar;
            m_progress_bar = NULL;
        }

        // Clean up

        delete filteredIndex;
    }

    template <typename T>
    void
    ScaleSpaceFilter<T>::apply_parallellized_on_dimension(FeatureSpace<T> *fs,
            ArrayIndex<T> *originalIndex,
            ArrayIndex<T> *filteredPoints,
            size_t fixedDimension)
    {
        using namespace std;

        const CoordinateSystem<T> *cs = fs->coordinate_system;
        const vector<size_t> dim_sizes = cs->get_dimension_sizes();
        const size_t numDims = cs->rank();
        const size_t value_rank = fs->value_rank();

        // Collate the dimension sizes for the linear mapping
        // and the dimension indexes

        std::vector<size_t> iter_dimensions;

        for (size_t i = 0; i < dim_sizes.size(); i++)
            if (i != fixedDimension)
                iter_dimensions.push_back(dim_sizes[i]);

        LinearIndexMapping mapping(iter_dimensions);
        size_t N = mapping.size();

#if WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t n = 0; n < N; n++) {
            // re-assemble the starting point for the iteration from
            // the mapping's components

            vector<int> truncatedGridPoint = mapping.linear_to_grid(n);
            vector<int> startPoint(cs->rank());

            size_t runningIndex = 0;
            for (size_t i = 0; i < cs->rank(); i++)
                if (i != fixedDimension)
                    startPoint[i] = truncatedGridPoint[runningIndex++];

            size_t fixedDimSize = dim_sizes[fixedDimension];

            // Now iterate over the remaining dimension
            // TODO: consider to nest this depending on 
            // problem size

            for (size_t k = 0; k < dim_sizes[fixedDimension]; k++) {
                if (this->show_progress()) {
#if WITH_OPENMP
#pragma omp critical
#endif
                    m_progress_bar->operator++();
                }

                // construct position

                vector<int> gridpoint(startPoint);
                gridpoint[fixedDimension] = k;

#if SCALE_SPACE_SKIPS_NON_ORIGINAL_POINTS
                if (fs->off_limits()->get(gridpoint))
                    continue;
#endif  

                ScaleSpaceKernel<T> g = this->m_kernels[fixedDimension];

                // Find the boundaries. Take care not to step
                // outside the bounds of the array

                int width = g.values().size() - 1;
                int gpIndex = (int) gridpoint[fixedDimension];
                int minIndex = (gpIndex - width >= 0) ? (gpIndex - width) : 0;
                int maxIndex = ((gpIndex + width) < (fixedDimSize - 1)) ? (gpIndex + width) : (fixedDimSize - 1);

                // Convolute in 1D around the given point with
                // the mask size determined by the kernel
                // Run the convolution for each feature variable

                typename CoordinateSystem<T>::GridPoint gridIter = gridpoint;

                vector<T> sum(value_rank, 0.0);

                size_t sumCount = 0;

                for (int i = minIndex; i < maxIndex; i++) {
                    gridIter[fixedDimension] = i;

                    // Again, make sure no points originally marked as
                    // off limits are used

#if SCALE_SPACE_SKIPS_NON_ORIGINAL_POINTS
                    if (fs->off_limits()->get(gridpoint))
                        continue;
#endif  

                    // get the point at the iterated position

                    Point<T> *pIter = originalIndex->get(gridIter);

                    if (pIter == NULL)
                        continue;

                    // apply the pre-sampled gaussian and sum up
                    // in each variable in the feature-space

                    size_t d = (i <= k) ? (k - i) : (i - k);

                    for (int varIndex = 0; varIndex < value_rank; varIndex++) {
                        T value = pIter->values[cs->rank() + varIndex];
                        sum[varIndex] += g.value(d) * value;
                    }

                    sumCount++;
                }

                // No muss, no fuss

                if (sumCount == 0)
                    continue;

                // Fuss! Fetch the point from the array index

                Point<T> *p = filteredPoints->get(gridpoint);

                // If no point existed, decide if we need to create one

                if (p == NULL) {
                    // Create a new point with default values
                    // and insert into array index

                    typename CoordinateSystem<T>::Coordinate coordinate = cs->newCoordinate();
                    cs->lookup(gridpoint, coordinate);
                    vector<T> values = coordinate;
                    values.resize(fs->rank(), 0.0);

                    p = PointFactory<T>::get_instance()->create(gridpoint, coordinate, values);

                    // Did this exist in the original index?

                    Point<T> *op = originalIndex->get(gridpoint);
                    p->isOriginalPoint = ((op == NULL) ? false : op->isOriginalPoint);

                    // Since we just created this point, there
                    // is no need to copy it again when adding
                    // it to the array index

#if WITH_OPENMP
#pragma omp critical
                    {
#endif
                    filteredPoints->set(gridpoint, p, false);
#if WITH_OPENMP
                    }
#endif
                }

                // If we have a point after all that, update it with the
                // filtered value

                if (p != NULL) {
                    // copy values and track limits

                    for (int varIndex = 0; varIndex < value_rank; varIndex++) {
                        p->values[numDims + varIndex] = sum[varIndex];

                        if (sum[varIndex] < m_min[varIndex]) {
#if WITH_OPENMP
#pragma omp critical
#endif
                            m_min[varIndex] = sum[varIndex];
                        }

                        if (sum[varIndex] > m_max[varIndex]) {
#if WITH_OPENMP
#pragma omp critical
#endif
                            m_max[varIndex] = sum[varIndex];
                        }
                    }
                }
            }
        }
    }

    template <typename T>
    void
    ScaleSpaceFilter<T>::apply_parallellized(FeatureSpace<T> *fs)
    {
        LinearIndexMapping mapping(fs->coordinate_system->get_dimension_sizes());
        size_t N = mapping.size();
        const size_t value_rank = fs->value_rank();

        using namespace std;

        const CoordinateSystem<T> *cs = fs->coordinate_system;

        // index the original

        if (this->show_progress()) {
            cout << endl << "Constructing array indexes ...";
            start_timer();
        }

        ArrayIndex<T> *originalIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), fs->points, true);
        ArrayIndex<T> *filteredIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), false);

        if (this->show_progress()) {
            cout << "done. (" << stop_timer() << "s)" << endl;
        }

        if (this->show_progress()) {
            cout << endl << "Applying scale filter t=" << m_scale << " decay=" << m_decay << " ... " << endl;

            long numPoints = 1;
            for (size_t i = 0; i < fs->coordinate_system->rank(); i++) {
                numPoints *= fs->coordinate_system->dimensions()[i].getSize();
            }
            m_progress_bar = new boost::progress_display(fs->spatial_rank() * numPoints);
            start_timer();
        }

        // initialize min/max and re-set counts
        for (size_t varIndex = 0; varIndex < value_rank; varIndex++) {
            m_min[varIndex] = std::numeric_limits<T>::max();
            m_max[varIndex] = std::numeric_limits<T>::min();
        }

        //
        // Apply dimension by dimension (exploiting separability)
        //

        for (size_t dimIndex = 0; dimIndex < fs->spatial_rank(); dimIndex++) {
            apply_parallellized_on_dimension(fs, originalIndex, filteredIndex, dimIndex);

            delete originalIndex;

            if (dimIndex < (fs->spatial_rank() - 1)) {
                originalIndex = filteredIndex;

                filteredIndex = new ArrayIndex<T>(cs->get_dimension_sizes(), false);
            }
        }

        // replace the points in the original with the filtered
        // array index results

        filteredIndex->replace_points(fs->points);

        size_t originalPoints = 0;
        for (size_t i = 0; i < fs->points.size(); i++) {
            if (fs->points[i]->isOriginalPoint) originalPoints++;
        }

        if (this->show_progress()) {
            cout << "done. (" << stop_timer() << "s)" << endl;
            cout << "Filtered featurespace contains " << fs->size() << " points (" << originalPoints << " original points, "
                    << "(" << (fs->size() - originalPoints) << " new points))" << endl;
            delete m_progress_bar;
            m_progress_bar = NULL;
        }

        // Clean up

        delete filteredIndex;
    }

    template <typename T>
    void
    ScaleSpaceFilter<T>::apply(FeatureSpace<T> *fs)
    {
        this->m_unfiltered_min = fs->min();
        this->m_unfiltered_max = fs->max();
        //this->applyWithArrayIndex(fs);

        this->apply_parallellized(fs);
    }

#pragma mark -
#pragma mark Range handling

    template <typename T>
    map<size_t, T>
    ScaleSpaceFilter<T>::getRangeFactors()
    {
        map<size_t, T> factors;
        typename map<size_t, T>::iterator mi;

        for (mi = m_min.begin(); mi != m_min.end(); mi++) {
            size_t i = mi->first;
            factors[i] = (m_max[i] - m_min[i]) / (m_unfiltered_max[i] - m_unfiltered_min[i]);
        }

        return factors;
    }

    template <typename T>
    const map<size_t, T> &
    ScaleSpaceFilter<T>::get_filtered_min()
    {
        return m_min;
    }

    template <typename T>
    const map<size_t, T> &
    ScaleSpaceFilter<T>::get_filtered_max()
    {
        return m_max;
    }
}

#endif
