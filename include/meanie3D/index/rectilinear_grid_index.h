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

#ifndef M3D_RECTILINEAR_GRID_INDEX_H
#define M3D_RECTILINEAR_GRID_INDEX_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/index.h>

namespace m3D {

    /** Implementation of index using ArrayIndex and searching by grid points 
     * rather than using a KD-Tree. Condition: rectilinear and (approximately)
     * uniform grid (in each dimension)
     */
    template <typename T>
    class RectilinearGridIndex : public PointIndex<T>
    {
        friend class PointIndex<T>;

    private:

#pragma mark -
#pragma mark Member variables

        ArrayIndex<T> *m_index;

    protected:

#pragma mark
#pragma mark Protected Constructor/Destructor

        inline
        RectilinearGridIndex(typename Point<T>::list *points, size_t dimension) : PointIndex<T>(points, dimension), m_index(NULL)
        {
        };

        inline
        RectilinearGridIndex(typename Point<T>::list *points, const vector<size_t> &indexes) : PointIndex<T>(points, indexes), m_index(NULL)
        {
        };

        inline
        RectilinearGridIndex(FeatureSpace<T> *fs) : PointIndex<T>(fs), m_index(NULL)
        {
        };

        inline
        RectilinearGridIndex(FeatureSpace<T> *fs, const vector<netCDF::NcVar> &index_variables) : PointIndex<T>(fs, index_variables), m_index(NULL)
        {
        };

#pragma mark -
#pragma mark Destructor

    public:

        ~RectilinearGridIndex()
        {
            if (m_index != NULL) {
                delete m_index;

                m_index = NULL;
            }
        };

#pragma mark -
#pragma mark Copy Operator

        /** Copy operator
         */
        RectilinearGridIndex<T> operator=(const RectilinearGridIndex<T> &other)
        {
            return RectilinearGridIndex<T>(other);
        }

        void
        build_index(const vector<T> &ranges)
        {
            // Dispose of previous indexes but inform user

            if (m_index != NULL) {
                cerr << "WARNING:called RectilinearGridIndex::build_index twice. Disposing of previous array index" << endl;

                delete m_index;

                m_index = NULL;
            }

            vector<size_t> dimensions = this->m_fs->coordinate_system->get_dimension_sizes();

            m_index = new ArrayIndex<T>(dimensions, this->m_fs->points, false);
        }

    public:

#pragma mark -
#pragma mark Overwritten Public Methods

        void
        remove_point(typename Point<T>::ptr p)
        {
            m_index->set(p->gridpoint, NULL);
        }

        void
        add_point(typename Point<T>::ptr p)
        {
            m_index->set(p->gridpoint, p);
        }

        typename Point<T>::list *
        search(const vector<T> &x, const SearchParameters *params, vector<T> *distances = NULL)
        {
            vector<T> h;

            if (params->search_type() == SearchTypeRange) {
                RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;

                h = p->bandwidth;
            } else {
                throw "RectilinearGridIndex does not support knn search at this time";
            }

            // spatial realm

            typename CoordinateSystem<T>::GridPoint gp = this->m_fs->coordinate_system->newGridPoint();

            vector<T> x_spatial = this->m_fs->spatial_component(x);

            typename Point<T>::list *result = NULL;

            try {
                this->m_fs->coordinate_system->reverse_lookup(x_spatial, gp);

                result = this->search(x, gp, h, distances);
            } catch (std::out_of_range& e) {
                cerr << "ERROR:reverse coordinate transformation failed for coordinate=" << x_spatial << endl;
            }

            return result;
        }

#pragma mark -
#pragma mark Public Methods (new)

        typename Point<T>::list *
        search(const vector<T> &x,
                const typename CoordinateSystem<T>::GridPoint &gp,
                const vector<T> &h,
                vector<T> *distances = NULL)
        {
            if (m_index == NULL) {
                this->build_index(h);
            }

            // Calculate the grid points

            typename Point<T>::list *result = new typename Point<T>::list();

            vector<size_t> lower_index_bounds(h.size(), 0);

            vector<size_t> upper_index_bounds(h.size(), 0);

            const CoordinateSystem<T> *cs = this->m_fs->coordinate_system;

            vector<T> resolution = cs->resolution();

            for (size_t i = 0; i < resolution.size(); i++) {
                size_t num_gridpoints = h[i] / resolution[i];

                // bound against underrun

                lower_index_bounds[i] = ((gp[i] - num_gridpoints) <= 0)
                        ? 0
                        : lower_index_bounds[i] = gp[i] - num_gridpoints;

                // bound against overrun

                size_t max_index = cs->dimensions()[i].getSize() - 1;

                upper_index_bounds[i] = ((gp[i] + num_gridpoints) >= max_index)
                        ? max_index
                        : upper_index_bounds[i] = gp[i] + num_gridpoints;
            }

            // Start the recursion

            typename CoordinateSystem<T>::GridPoint gridpoint = cs->newGridPoint();

            search_recursive(0, x, gridpoint, h, lower_index_bounds, upper_index_bounds, result, distances);

            if (PointIndex<T>::write_index_searches) {
                this->write_search(x, h, result);
            }

            return result;
        }

#pragma mark
#pragma mark Protected specific methods

    protected:

        void
        search_recursive(size_t dim_index,
                const vector<T> &x,
                typename CoordinateSystem<T>::GridPoint &gridpoint,
                const vector<T> &h,
                const vector<size_t> &lower_index_bounds,
                const vector<size_t> &upper_index_bounds,
                typename Point<T>::list *result,
                vector<T> *distances = NULL)
        {
            if (dim_index < (gridpoint.size() - 1)) {
                for (size_t index = lower_index_bounds[dim_index]; index <= upper_index_bounds[dim_index]; index++) {
                    gridpoint[dim_index] = index;

                    search_recursive(dim_index + 1, x, gridpoint, h, lower_index_bounds, upper_index_bounds, result, distances);
                }
            } else {
                for (size_t index = lower_index_bounds[dim_index]; index <= upper_index_bounds[dim_index]; index++) {
                    gridpoint[dim_index] = index;

                    typename Point<T>::ptr p = m_index->get(gridpoint);

                    if (p != NULL) {
                        // value realm will have to be calculated with a
                        // simple euclidean metric for now

                        bool within_range = true;

                        vector<T> dists(x.size(), 0.0);

                        for (size_t hi = gridpoint.size(); (hi < h.size()) && within_range; hi++) {
                            T dist = boost::numeric_cast<T>(abs(x[hi] - p->values[hi]));

                            if (distances != NULL) {
                                dists[hi] = dist;
                            }

                            within_range = (dist <= h[hi]);
                        }

                        if (within_range) {
                            result->push_back(p);

                            if (distances != NULL) {
                                // fill out the spatial realm part of the dists vector

                                for (size_t i = 0; i < gridpoint.size(); i++) {
                                    dists[i] = boost::numeric_cast<T>(abs(x[i] - p->values[i]));
                                }

                                // distances is a component stretching across the whole realm
                                // in what way this still makes sense, must be questioned. In
                                // a future implementation the distance should be given for
                                // each component individually?

                                distances->push_back(utils::vectors::vector_norm<T>(dists));
                            }
                        }
                    }
                }
            }
        }
    };
}

#endif
