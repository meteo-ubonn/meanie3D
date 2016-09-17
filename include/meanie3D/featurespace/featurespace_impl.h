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


#ifndef M3D_FEATURESPACE_IMPL_H
#define M3D_FEATURESPACE_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>

#include <algorithm>
#include <limits>
#include <exception>

#include <boost/filesystem.hpp>

#include "featurespace.h"

namespace m3D {

    template <typename T>
    map<int, double> *FeatureSpace<T>::NO_THRESHOLDS = NULL;

    template <typename T>
    class MeanshiftOperation;

#pragma mark -
#pragma mark Constructors

    template <typename T>
    FeatureSpace<T>::FeatureSpace(const CoordinateSystem<T> *coordinate_system,
            const DataStore<T> *data_store,
            const map<int, double> &lower_thresholds,
            const map<int, double> &upper_thresholds,
            const map<int, double> &replacement_values,
            const bool& show_progress)
    : m_data_store(data_store)
    , m_progress_bar(NULL)
    , m_lower_thresholds(lower_thresholds)
    , m_upper_thresholds(upper_thresholds)
    , m_replacement_values(replacement_values)
    , m_off_limits(NULL)
    , coordinate_system(coordinate_system)
    {
        dimension = coordinate_system->rank() + data_store->rank();

        // construct feature space
        this->construct_featurespace(show_progress);
    }

    /** Destructor
     */
    template <typename T>
    FeatureSpace<T>::~FeatureSpace()
    {
        if (m_off_limits != NULL) {
            delete m_off_limits;
        }
    }

#pragma mark -
#pragma mark Copy code

    template <typename T>
    FeatureSpace<T>::FeatureSpace(const FeatureSpace<T>& other, bool with_points)
    : m_progress_bar(NULL)
    , m_data_store(other.data_store())
    , coordinate_system(other.coordinate_system)
    , m_lower_thresholds(other.m_lower_thresholds)
    , m_upper_thresholds(other.m_upper_thresholds)
    , m_off_limits(other.m_off_limits)
    {
        if (with_points) {
            // Copy points

            typename Point<T>::list::const_iterator pi;

            for (pi = other.points.begin(); pi != other.points.end(); pi++) {
                typename Point<T>::ptr p = *pi;
                typename Point<T>::ptr copy = PointFactory<T>::get_instance()->copy(p);
                points.push_back(copy);
            }
        }
    }

    template <typename T>
    FeatureSpace<T>::FeatureSpace(const FeatureSpace<T>* other,
            bool with_points)
    : m_progress_bar(NULL)
    , m_data_store(other->data_store())
    , coordinate_system(other->coordinate_system)
    , m_lower_thresholds(other->m_lower_thresholds)
    , m_upper_thresholds(other->m_upper_thresholds)
    , m_off_limits(other->m_off_limits)
    {
        if (with_points) {
            typename Point<T>::list::const_iterator pi;

            for (pi = other->points.begin(); pi != other->points.end(); pi++) {
                typename Point<T>::ptr p = *pi;
                typename Point<T>::ptr copy = PointFactory<T>::get_instance()->copy(p);
                points.push_back(copy);
            }
        }
    }

    // Assignment Operator =

    template <typename T>
    FeatureSpace<T>
    FeatureSpace<T>::operator=(const FeatureSpace<T>& other)
    {
        return FeatureSpace<T>(other);
    }

#pragma mark -
#pragma mark Building feature space maps

    template <typename T>
    void FeatureSpace<T>::construct_featurespace(bool show_progress)
    {
        if (show_progress) {
            cout << endl << "Constructing feature space ... ";
        }

        using boost::progress_display;

        // Iterate over variables

        if (show_progress) {
            size_t number_of_points = 1;

            for (size_t di = 0; di < coordinate_system->rank(); di++) {
                NcDim dim = coordinate_system->dimensions()[di];
                number_of_points *= dim.getSize();
            }

            m_progress_bar = new progress_display(number_of_points);
        }

        // initialize min/max

        for (size_t vi = 0; vi < data_store()->rank(); vi++) {
            m_min[vi] = std::numeric_limits<T>::max();
            m_max[vi] = std::numeric_limits<T>::min();
        }

        this->build();

        if (show_progress) {
            cout << " done. (" << size() << " points in " << stop_timer() << " seconds)" << endl;
        }

        if (show_progress) {
            delete m_progress_bar;
        }
    }

    template <typename T>
    void FeatureSpace<T>::build()
    {
        m_off_limits = new MultiArrayBlitz<bool>(this->coordinate_system->get_dimension_sizes(), false);

        size_t size = this->m_data_store->size();

        LinearIndexMapping mapping(m_data_store->coordinate_system()->get_dimension_sizes());

#if WITH_OPENMP
#pragma omp parallel for schedule(dynamic,10) 
#endif
        for (size_t linear_index = 0; linear_index < size; linear_index++) {
            // Get the variables together and construct the cartesian coordinate
            if (m_progress_bar != NULL) {
#if WITH_OPENMP
#pragma omp critical
#endif
                m_progress_bar->operator++();
            }

            vector<int> gridpoint = mapping.linear_to_grid(linear_index);

            typename CoordinateSystem<T>::Coordinate coordinate(gridpoint.size());

            coordinate_system->lookup(gridpoint, coordinate);

            // Iterate over the variables

            bool isPointValid = true;
            bool isPointInRange = true;

            // start the entry by copying the dimension variables

            vector<T> values = coordinate;

            // iterate over the variables

            for (size_t var_index = 0; var_index < this->m_data_store->rank() && isPointValid; var_index++) {
                typename std::map<int, double>::const_iterator replacement;
                replacement = this->m_replacement_values.find(var_index);

                // is this contribution valid?

                T value = data_store()->get(var_index, gridpoint, isPointInRange, isPointValid);

                if (!isPointValid) {
                    // Reading routine marked this point 'off limits'
#if WITH_OPENMP
#pragma omp critical
#endif
                    this->m_off_limits->set(gridpoint, true);
                }

                if (isPointValid && isPointInRange) {
                    values.push_back(value);

                    // apply upper/lower thresholding to the value, if asked

                    map<int, double>::const_iterator fi;
                    fi = m_lower_thresholds.find(var_index);
                    if (fi != m_lower_thresholds.end()) {
                        if (value < fi->second) {
                            if (replacement != this->m_replacement_values.end()) {
                                value = replacement->second;
                            } else {
                                isPointValid = false;
                            }
                        }
                    }

                    fi = m_upper_thresholds.find(var_index);
                    if (fi != m_upper_thresholds.end()) {
                        if (value > fi->second) {
                            if (replacement != this->m_replacement_values.end()) {
                                value = replacement->second;
                            } else {
                                isPointValid = false;
                            }
                        }
                    }
                } else if (replacement != this->m_replacement_values.end()) {
                    // check for replacement value and use after all
                    value = replacement->second;
                    isPointValid = true;
                }
            }

            // if the point is still valid (all variables measured up to criteria)
            // add it to the feature-space

            if (isPointValid) {
                typename Point<T>::ptr p = NULL;

#if WITH_OPENMP
#pragma omp critical
#endif
                {
                    p = PointFactory<T>::get_instance()->create(gridpoint, coordinate, values);
                    p->isOriginalPoint = true;
                    this->points.push_back(p);
                }

                for (size_t vi = 0; vi < data_store()->rank(); vi++) {
                    if (p->values[coordinate.size() + vi] < m_min[vi]) {
#if WITH_OPENMP
#pragma omp critical
#endif
                        m_min[vi] = p->values[coordinate.size() + vi];
                    }

                    if (p->values[coordinate.size() + vi] > m_max[vi]) {
#if WITH_OPENMP
#pragma omp critical
#endif
                        m_max[vi] = p->values[coordinate.size() + vi];
                    }
                }
            }
        }
    }

#pragma mark -
#pragma mark Other

    template <typename T>
    typename CoordinateSystem<T>::Coordinate
    FeatureSpace<T>::spatial_component(const vector<T> &value) const
    {
        assert(value.size() >= this->coordinate_system->rank());

        typename CoordinateSystem<T>::Coordinate coordinate(&value[0], &value[ this->coordinate_system->rank() ]);

        return coordinate;
    }

    template <typename T>
    T
    FeatureSpace<T>::get_spatial_component_at(typename Point<T>::ptr p, size_t index) const
    {
        assert(index < this->coordinate_system->rank());
    }

    template <typename T>
    T
    FeatureSpace<T>::get_value_component_at(typename Point<T>::ptr p, size_t index) const
    {
        assert(index < this->data_store()->rank());

    }

    template <typename T>
    void
    FeatureSpace<T>::round_to_resolution(vector<T> &point,
            const vector<T> &resolution) const
    {
        assert(point.size() <= resolution.size());

        for (size_t i = 0; i < point.size(); i++) {
#if GRID_ROUNDING_METHOD_FLOOR
            T multiplier = floor(point[i] / resolution[i]);
#elif GRID_ROUNDING_METHOD_ROUND
            T multiplier = round(point[i] / resolution[i]);
#elif GRID_ROUNDING_METHOD_CEIL
            T multiplier = ceil(point[i] / resolution[i]);
#elif GRID_ROUNDING_METHOD_RINT
            T multiplier = rint(point[i] / resolution[i]);
#elif GRID_ROUNDING_METHOD_NONE
            T multiplier = point[i] / resolution[i];
#else
            T multiplier = ceil(point[i] / resolution[i]);
#endif
            point[i] = multiplier * resolution[i];
        }
    }

    template <typename T>
    void
    FeatureSpace<T>::print() const
    {
        using namespace std;

        cout << "\t#points = " << size() << endl;

        for (size_t i = 0; i < points.size(); i++) {
            typename Point<T>::ptr p = points[i];

            cout << "\t\t#" << i << "\t" << p->values << endl;
        }
    }

    template <typename T>
    void
    FeatureSpace<T>::clear()
    {
        for (size_t i = 0; i < points.size(); i++) {
            typename Point<T>::ptr p = points[i];
            points[i] = NULL;
            delete p;
        }
        points.clear();
    }

    template <typename T>
    size_t
    FeatureSpace<T>::count_original_points() const
    {
        size_t originalPoints = 0;
        for (size_t i = 0; i < this->points.size(); i++) {
            if (this->points[i]->isOriginalPoint) originalPoints++;
        }
        return originalPoints;
    }

    template <typename T>
    void
    FeatureSpace<T>::sanity_check()
    {
        size_t originalPoints = 0;
        for (size_t i = 0; i < this->points.size(); i++) {
            typename Point<T>::ptr p = this->points[i];

            if (p->gridpoint.size() != this->spatial_rank()
                    || p->coordinate.size() != this->spatial_rank()
                    || p->values.size() != this->rank()) {
                cerr << "Point #" << i << " is insane:"
                        << " gridpoint.size()=" << p->gridpoint.size()
                        << " coordinate.size()=" << p->coordinate.size()
                        << " values.size()=" << p->values.size()
                        << endl;

                throw std::range_error("insane point");
            }
        }
    }

}

#endif
