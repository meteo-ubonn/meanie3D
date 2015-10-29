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

#ifndef M3D_LINEARINDEX_H
#define M3D_LINEARINDEX_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/index.h>
#include <meanie3D/featurespace/point.h>

#include <iostream>
#include <algorithm>

namespace m3D {

    /** Implementation of FeatureSpace which simply searches the feature-space vector
     * brute-force style when sampling around points.
     * 
     * TODO: this implementation is NOT UP TO DATE
     */
    template <typename T>
    class LinearIndex : public PointIndex<T>
    {
        friend class PointIndex<T>;

    protected:

#pragma mark -
#pragma mark Constructor/Destructor

        inline
        LinearIndex(typename Point<T>::list *points, size_t dimension) : PointIndex<T>(points, dimension)
        {
        };

        inline
        LinearIndex(typename Point<T>::list *points, const vector<size_t> &indexes) : PointIndex<T>(points, indexes)
        {
        };

        inline
        LinearIndex(FeatureSpace<T> *fs) : PointIndex<T>(fs)
        {
        };

        inline
        LinearIndex(FeatureSpace<T> *fs, const vector<netCDF::NcVar> &index_variables) : PointIndex<T>(fs, index_variables)
        {
        };

        inline
        LinearIndex(const LinearIndex<T> &o) : PointIndex<T>(o)
        {
        };

#pragma mark -
#pragma mark Overwritten Protected Methods

        void
        build_index(const vector<T> &ranges)
        {
            // nothing to do. Index is the feature space itself

            // TODO: construct feature space from index variables
        };

    public:

        inline
        ~LinearIndex()
        {
        }

#pragma mark -
#pragma mark Copy Operator

        /** Copy operator
         */
        LinearIndex<T> operator=(const LinearIndex<T> &other)
        {
            return LinearIndex<T>(other);
        }

#pragma mark -
#pragma mark Overwritten Public Methods

        typename Point<T>::list *
        search(const vector<T> &x, const SearchParameters *params, vector<T> *distances = NULL)
        {
            using std::cerr;
            using std::endl;

            typename Point<T>::list *result = new typename Point<T>::list();

            // Test if a point is within the given ellipsoid

            // x1^2/a1^2 + .... + xn^2/an^2 <= 1

            // pre-calculate the denominators

            if (params->search_type() == SearchTypeKNN) {
                //                KNNSearchParams *p = (KNNSearchParams *)  &params;

                cerr << "FATAL:KNN is not supported yet in LinearIndex" << endl;
                exit(EXIT_FAILURE);
            }

            RangeSearchParams<T> *p = (RangeSearchParams<T> *) & params;

            vector<T> coefficients(p->bandwidth);

            for (size_t index = 0; index < p->bandwidth.size(); index++) {
                coefficients[index] = 1.0 / (coefficients[index] * coefficients[index]);
            }

            typename Point<T>::list::const_iterator it = this->m_fs->points.begin();

            while (it != this->m_fs->points.end()) {
                typename Point<T>::ptr p = *it;

                // r = x1^2/a1^2 + .... + xn^2/an^2

                float r = 0.0;

                vector<T> coordinate = p->coordinate;

                for (size_t index = 0; index < coordinate.size(); index++) {
                    float dist = coordinate[index] - x[index];

                    r += (dist * dist * coefficients[index]);
                }

                if (r <= 1.0) {
                    result->push_back(p);
                }

                it++;
            }

            return result;
        };

        void
        add_point(typename Point<T>::ptr p)
        {
            this->m_fs->points.push_back(p);
        }

        void
        remove_point(typename Point<T>::ptr p)
        {
            typename Point<T>::list::iterator f = find(this->m_fs->points.begin(), this->m_fs->points.end(), p);

            if (f != this->m_fs->points.end()) {
                this->m_fs->points.erase(f);
            }
        }
    };
}

#endif
