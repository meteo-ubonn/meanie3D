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

#ifndef M3D_FS_INDEXED_FLANN_H
#define M3D_FS_INDEXED_FLANN_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/index.h>

#include <flann/flann.hpp>
#include <flann/io/hdf5.h>

namespace m3D {

    using flann::Matrix;
    using flann::Index;
    using flann::L2;

    /** Implementation of FeatureSpace which simply searches the feature-space vector
     * brute-force style when sampling around points.
     */
    template <typename T>
    class FLANNIndex : public WhiteningIndex<T>
    {
        friend class PointIndex<T>;

    private:

#pragma mark -
#pragma mark Member variables

        Index< L2<T> > *m_index;
        Matrix<T> m_dataset;

    protected:

#pragma mark
#pragma mark Protected Constructor/Destructor

        inline
        FLANNIndex(typename Point<T>::list *points, size_t dimension) : WhiteningIndex<T>(points, dimension), m_index(NULL)
        {
        };

        inline
        FLANNIndex(typename Point<T>::list *points, const vector<size_t> &indexes) : WhiteningIndex<T>(points, indexes), m_index(NULL)
        {
        };

        inline
        FLANNIndex(FeatureSpace<T> *fs) : WhiteningIndex<T>(fs), m_index(NULL)
        {
        };

        inline
        FLANNIndex(FeatureSpace<T> *fs, const vector<netCDF::NcVar> &index_variables) : WhiteningIndex<T>(fs, index_variables), m_index(NULL)
        {
        };

        inline
        FLANNIndex(const FLANNIndex<T> &o) : WhiteningIndex<T>(o), m_dataset(dynamic_cast<FLANNIndex> (o).dataset())
        {
            build_index_from_dataset();
        };

#pragma mark -
#pragma mark Destructor

    public:

        ~FLANNIndex()
        {
            if (m_index != NULL) {
                delete m_index;

                m_index = NULL;
            }

            if (m_dataset.ptr()) {
                delete[] m_dataset.ptr();
            }
        };

#pragma mark -
#pragma mark Copy Operator

        /** Copy operator
         */
        FLANNIndex<T> operator=(const FLANNIndex<T> &other)
        {
            return FLANNIndex<T>(other);
        }

#pragma mark -
#pragma mark Overwritten Protected Methods

    protected:

        void
        build_index(const vector<T> &ranges)
        {
            this->transform_featurespace(ranges);

            if (m_index != NULL) {
                delete m_index;

                m_index = NULL;
            }

            // Re-package data for FLANN
            construct_dataset();

            build_index_from_dataset();

        }

    public:

#pragma mark -
#pragma mark Overwritten Public Methods

        void
        remove_point(typename Point<T>::ptr p)
        {
            throw "FLANN does not support removing individual points";
        }

        void
        add_point(typename Point<T>::ptr p)
        {
            // construct the coordinate for indexing

            vector<T> coordinate(this->size());

            for (size_t i = 0; i<this->dimension(); i++) {
                size_t j = this->m_index_variable_indexes[i];

                coordinate[i] = p->coordinate[j];
            }

            // transform it into the whitened space

            vector<T> t = this->transform_vector(coordinate);

            // add it to the index

            T data[t.size()][1];

            for (size_t i = 0; i<this->dimension(); i++) {
                data[i][0] = t[i];
            }

            Matrix<T> dataset(&data[0][0], 1, t.size());

            m_index->addPoints(dataset);

            // add to the points list

            this->m_points->push_back(p);
        }

        typename Point<T>::list *
        search(const vector<T> &x, const SearchParameters *params, vector<T> *distances = NULL)
        {
            // Check if the index needs re-building

            vector<T> h;

            if (params->search_type() == SearchTypeRange) {
                RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;

                h = p->bandwidth;
            } else {
                h = vector<T>(x.size(), 1.0);
            }

            if (this->white_range != h || m_index == NULL) {
                build_index(h);
            }

            // FLANN Search Parameters

            flann::SearchParams flann_params(32, 0, false);
            flann_params.use_heap = flann::FLANN_True;
            flann_params.matrices_in_gpu_ram = flann::FLANN_True;

            // Build query

            vector<T> x_t = (params->search_type() == SearchTypeRange) ? this->transform_vector(x) : x;

            T query_point[x_t.size()];

            for (size_t i = 0; i < x_t.size(); i++) {
                query_point[i] = x_t[i];
            }

            flann::Matrix<T> query = flann::Matrix<T>(&query_point[0], 1, this->dimension());

#if DEBUG_INDEX_SEARCHES
            double t = utils::stop_timer();
            std::cout << std::endl << "\tsearching around x = " << x
                    << " ->  " << x_t << " ...";
            utils::start_timer();
#endif
            // Go get it!

            vector< vector<int> > indices;

            vector<vector<T> > dists;

            if (params->search_type() == SearchTypeKNN) {
                const KNNSearchParams<T> *p = dynamic_cast<const KNNSearchParams<T> *> (params);

                flann_params.sorted = flann::FLANN_True;

                m_index->knnSearch(query, indices, dists, p->k, flann_params);
            } else {
                m_index->radiusSearch(query, indices, dists, this->white_radius, flann_params);
            }

#if DEBUG_INDEX_SEARCHES
            std::cout << " with radius = " << this->white_radius;
            std::cout << " found " << indices[0].size() << " points. (" << t << "s)" << std::endl;
#endif
            // re-wrap results

            typename Point<T>::list * result = new typename Point<T>::list();

            for (size_t row = 0; row < indices.size(); row++) {
                vector <int> _indices = indices[row];
                for (size_t col = 0; col < _indices.size(); col++) {
                    typename Point<T>::ptr p = this->m_points->at(_indices[col]);
                    result->push_back(p);
                }
            }

            if (distances) {
                for (size_t row = 0; row < dists.size(); row++) {
                    vector <T> _dists = dists[row];
                    for (size_t col = 0; col < _dists.size(); col++) {
                        distances->push_back((T) _dists[ col ]);
                    }
                }
            }

            if (PointIndex<T>::write_index_searches) {

                if (params->search_type() == SearchTypeRange) {
                    RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;
                    this->write_search(x, p->bandwidth, result);
                } else {
                    vector<T> white_ranges(x_t.size(), this->white_radius);
                    this->write_search(x_t, white_ranges, result);
                }
            }

            return result;
        }

#pragma mark
#pragma mark Protected specific methods

    protected:

        /** Protected accessor. Required for copy constructor
         * @return pointer to the index
         */
        Index< L2<T> > *index()
        {
            return m_index;
        };

        /** Protected accessor. Required for copy constructor
         * @return dataset
         */
        Matrix<T> dataset()
        {
            return m_dataset;
        };

        /** Create dataset from the feature-space, as prescribed
         * by the index variables
         */
        void construct_dataset()
        {
            // re-package data for FLANN

            size_t count = this->size() * this->dimension();

            T* data = (T *) malloc(count * sizeof (T));

            for (size_t row = 0; row < this->white_point_matrix.size1(); row++) {
                for (size_t col = 0; col < this->white_point_matrix.size2(); col++) {
                    data[ row * this->white_point_matrix.size2() + col ] = this->white_point_matrix(row, col);
                }
            }

            m_dataset = flann::Matrix<T>(data, this->size(), this->dimension());

#if DEBUG_INDEX
#if WRITE_INDEX
            static size_t index_count = 0;
            string fn = "flann_dataset_" + boost::lexical_cast<string>(index_count) + ".h5";
            std::cout << "Saving dataset to " << fn << " ... ";
            flann::save_to_file(m_dataset, fn, "whitened_fs");
            std::cout << "done." << std::endl;
#endif
#endif
        }

        void build_index_from_dataset()
        {
            // Set up Index options
            flann::IndexParams params;

            params = flann::KDTreeSingleIndexParams();

            m_index = new flann::Index< flann::L2<T> >(m_dataset, params);

            m_index->buildIndex();

#if DEBUG_INDEX
#if WRITE_INDEX
            static size_t index_count = 0;
            string fn = "flann_index_" + boost::lexical_cast<string>(index_count) + ".h5";
            std::cout << "Saving index to " << fn << " ... ";
            m_index->save(fn);
            std::cout << "done." << std::endl;
            index_count++;
#endif
#endif
        }
    };
}

#endif
