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

#ifndef M3D_REPLACEMENT_FILTER_IMPL_H
#define M3D_REPLACEMENT_FILTER_IMPL_H

#include <meanie3D/utils.h>
#include <meanie3D/index.h>
#include "replacement_filter.h"

#include <vector>

namespace m3D {

    using namespace std;

    template<typename T>
    ReplacementFilter<T>::ReplacementFilter(const ReplacementMode mode,
                                            const size_t variable_index,
                                            const std::vector<T> &bandwidth,
                                            const float percentage,
                                            bool show_progress)
            : FeatureSpaceFilter<T>(show_progress),
              m_replacement_mode(mode),
              m_variable_index(variable_index),
              m_bandwidth(bandwidth),
              m_percentage(percentage) {};

    template<typename T>
    void ReplacementFilter<T>::apply(FeatureSpace <T> *fs) {

        // Create a spatial index for the copied feature space
        PointIndex<T> *index = PointIndex<T>::create(fs->get_points(), fs->rank());
        SearchParameters *params = new RangeSearchParams<T>(m_bandwidth);
        size_t value_index = fs->spatial_rank() + m_variable_index;

        vector<T> filteredValues;
        filteredValues.resize(fs->size());

        typename Point<T>::list points;

#if WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
        for (size_t i = 0; i < fs->size(); i++) {

            // Get the values around the point (original index)
            T result = fs->points[i]->values[value_index];

            typename Point<T>::list *neighbours = NULL;

#if WITH_OPENMP
#pragma omp critical
            {
#endif
            neighbours = index->search(fs->points[i]->coordinate, params);
#if WITH_OPENMP
            }
#endif
            if (!(neighbours == NULL || neighbours->size() == 0)) {

                vector<T> values;
                for (size_t n = 0; n < neighbours->size(); n++) {
                    values.push_back(neighbours->at(n)->values[value_index]);
                }

                switch (m_replacement_mode) {

                    case ReplaceWithLowest:
                    case ReplaceWithMedian:
                        // sort the data in ascending order
                        std::sort(values.begin(), values.end());
                        break;

                    case ReplaceWithHighest:
                        // sort the data in descending order
                        std::sort(values.begin(), values.end(), std::greater<int>());
                        break;
                }

                // calculate the number of values that make up
                // the required percentage
                int num_values = round(values.size() * m_percentage);
                if (m_replacement_mode == ReplaceWithMedian) {
                    result = values[num_values / 2];
                } else {
                    // obtain the average of the last num_values values
                    T sum = 0.0;
                    for (int i = 0; i < num_values; i++)
                        sum += values[i];
                    result = sum / ((T) num_values);
                }
            }

            filteredValues[i] = result;
        }

        // Replace the values in the featurespace's points
        for (size_t i = 0; i < fs->size(); i++) {
            fs->points[i]->values[value_index] = filteredValues[i];
        }

        delete index;
        delete params;
    }
}

#endif
