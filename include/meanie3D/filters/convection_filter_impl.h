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

#ifndef M3D_REFLECTIVITYCONVECTIONFILTER_IMPL_H
#define M3D_REFLECTIVITYCONVECTIONFILTER_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/progress.hpp>

#include "convection_filter.h"

namespace m3D {

    template<typename T>
    ConvectionFilter<T>::ConvectionFilter(const vector <T> &bandwidth,
                                          const size_t index_of_z,
                                          const bool show_progress,
                                          const T z_convective,
                                          const T critical_delta_z,
                                          const T convective_radius_factor,
                                          const bool erase_non_convective)
            : FeatureSpaceFilter<T>(show_progress), m_bandwidth(bandwidth), m_index_of_z(index_of_z),
              m_convective_threshold(z_convective), m_critical_delta_z(critical_delta_z),
              m_convective_radius_factor(convective_radius_factor), m_erase_non_convective(erase_non_convective) {
    }

    template<typename T>
    ConvectionFilter<T>::~ConvectionFilter() {
    }

#pragma mark -
#pragma mark Abstract filter method

    template<typename T>
    void
    ConvectionFilter<T>::apply(FeatureSpace <T> *fs) {
        using namespace std;
        using namespace m3D::utils::vectors;

        boost::progress_display *progress_bar = NULL;

        if (this->show_progress()) {
            cout << endl << "Applying convection filter ...";

            progress_bar = new boost::progress_display(2 * fs->size());

            start_timer();
        }

        // Create a point index with spatial components only

        PointIndex<T> *index = PointIndex<T>::create(fs->get_points(), fs->rank());

        PointIndex<T> *convective_radius_index = PointIndex<T>::create(fs->get_points(), fs->rank());

        // params for index search

        RangeSearchParams<T> params(this->m_bandwidth);

        // convective 'radius'

        RangeSearchParams<T> convective_radius_params(m_convective_radius_factor *m_bandwidth);

        // Create a field to hold the convection mask

        MultiArrayBlitz<bool> convective_mask(fs->coordinate_system->get_dimension_sizes(), false);

        // Iterate over the feature-space to create the convective mask

        typename Point<T>::list::iterator pit;

        for (size_t k = 0; k < fs->points.size(); k++) {
            if (this->show_progress()) {
                progress_bar->operator++();
            }

            Point<T> *p = fs->points.at(k);

            typename Point<T>::list *sample = NULL;

            sample = index->search(p->coordinate, &params, NULL);

            // Check if the convective threshold is met

            bool is_convective = false;

            T z_p = p->values.at(m_index_of_z);

            if (z_p >= m_convective_threshold) {
                is_convective = true;
            } else {
                // Obtain the background reflectivity for this point

                T z_background = 0.0;

                for (size_t i = 0; i < sample->size(); i++) {
                    Point<T> *p = sample->at(i);

                    T z = p->values.at(m_index_of_z);

                    z_background += z;
                }

                // linear average

                z_background = z_background / boost::numeric_cast<T>(sample->size());

                // Now figure if this point classifies as 'convective' according to
                // the scheme

                // implicit: if z_background < 0
                // implicit: z_background < m_convective_threshold

                T deltaZ = 10.0;

                if (z_background >= 0 && z_background <= m_convective_threshold) {
                    deltaZ = 10.0 - (z_background * z_background) / 180.0;
                }

                if (deltaZ > m_critical_delta_z) {
                    is_convective = true;
                }
            }

            delete sample;

            if (is_convective) {
                // Mark all points within convective radius

                typename Point<T>::list *sample = NULL;

                sample = convective_radius_index->search(p->coordinate, &convective_radius_params, NULL);

                for (size_t i = 0; i < sample->size(); i++) {
                    Point<T> *sp = sample->at(i);

                    convective_mask.set(sp->gridpoint, true);
                }

                delete sample;
            }
        }

        // Don't need the index anymore

        delete index;
        delete convective_radius_index;

        // Now iterate over the feature-space again and set
        // reflectivity to zero for all non-convective points

        vector<Point<T> *> accepted;
        vector<Point<T> *> erased;

        for (size_t k = 0; k < fs->points.size(); k++) {
            if (this->show_progress()) {
                progress_bar->operator++();
            }

            Point<T> *p = fs->points.at(k);

            if (m_erase_non_convective) {
                if (convective_mask.get(p->gridpoint)) {
                    accepted.push_back(p);
                } else {
                    erased.push_back(p);
                }
            } else {
                p->values.at(m_index_of_z) = 0.0;
            }
        }

        if (m_erase_non_convective) {
            // replace feature-space points with those accepted
            // only and delete the rest

            fs->points = accepted;

            for (size_t i = 0; i < erased.size(); i++) {
                Point<T> *p = erased[i];
                delete p;
            }
        }

        if (this->show_progress()) {
            delete progress_bar;
            cout << "done. (" << stop_timer() << "s)" << std::endl;
        }
    }
}

#endif
