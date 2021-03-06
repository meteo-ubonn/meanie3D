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

#ifndef M3D_FEATURESPACE_KDTREE_H
#define M3D_FEATURESPACE_KDTREE_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/index.h>

#include <meanie3D/kdtree/kdtree.h>

namespace m3D {

    /** Implementation of FeatureSpace which simply searches the feature-space vector
     * brute-force style when sampling around points.
     */
    template<typename T>
    class KDTreeIndex : public WhiteningIndex<T>
    {
        friend class PointIndex<T>;

    private:

        // kdtree.c
        struct kdtree *m_index;

#pragma mark -
#pragma mark Member variables

    protected:

#pragma mark -
#pragma mark Constructor/Destructor

        inline
        KDTreeIndex(typename Point<T>::list *points, size_t dimension) : WhiteningIndex<T>(points, dimension),
                                                                         m_index(NULL) {
        };

        inline
        KDTreeIndex(typename Point<T>::list *points, const vector<size_t> &indexes) : WhiteningIndex<T>(points,
                                                                                                        indexes),
                                                                                      m_index(NULL) {
        };

        inline
        KDTreeIndex(FeatureSpace<T> *fs) : WhiteningIndex<T>(fs), m_index(NULL) {
        };

        inline
        KDTreeIndex(FeatureSpace<T> *fs, const vector<netCDF::NcVar> &index_variables) : WhiteningIndex<T>(fs,
                                                                                                           index_variables),
                                                                                         m_index(NULL) {
        };

        inline
        KDTreeIndex(const KDTreeIndex<T> &o) : WhiteningIndex<T>(o), m_index(dynamic_cast<KDTreeIndex> (o).kd_tree()) {
        };

#pragma mark -
#pragma mark Specific Protected Methods

        /** protected accessor for copy constructor 
         * @return kdtree 
         */
        struct kdtree *kd_tree() {
            return m_index;
        }

#pragma mark -
#pragma mark Overwritten Protected Methods

        void
        build_index(const vector<T> &ranges) {
            // TODO: build index from index variables

            this->transform_featurespace(ranges);

            if (m_index == NULL) {
                // allocate kd tree, if necessary
                m_index = kd_create(this->dimension());
            } else {
                // if there is a tree, clear it
                kd_clear(m_index);
            }

            typename Point<T>::list::iterator it;

            // Allocate 

            for (size_t row = 0; row < this->white_point_matrix.size1(); row++) {
                // re-package data for kd-tree

                double query_point[this->white_point_matrix.size2()];

                for (size_t col = 0; col < this->dimension(); col++) {
                    query_point[col] = this->white_point_matrix(row, col);
                }

                // Insert

                typename Point<T>::ptr p = this->m_points->at(row);

                kd_insert(m_index, query_point, static_cast<void *> (p));
            }
        }

    public:

#pragma mark -
#pragma mark Destructor

        /** Destructor 
         */
        inline
        ~KDTreeIndex() {
            if (m_index != NULL) {
                kd_free(m_index);
            }
        }

#pragma mark -
#pragma mark Copy Operator

        /** Copy operator 
         */
        KDTreeIndex<T> operator=(const KDTreeIndex<T> &other) {
            return KDTreeIndex<T>(other);
        }

#pragma mark -
#pragma mark Overwritten Public Methods   

        void
        remove_point(typename Point<T>::ptr p) {
            throw "KDTree does not support removing individual points";
        }

        void
        add_point(typename Point<T>::ptr p) {
            throw "KDTree does not support adding individual points";
        }

        typename Point<T>::list *
        search(const vector<T> &x, const SearchParameters *params, vector<T> *distances = NULL) {
            if (params->search_type() == SearchTypeKNN) {
                std::cerr << "FATAL:KNN not supported by KDTree yet" << std::endl;
                exit(EXIT_FAILURE);
            }

            // Check if the index needs re-building

            if (params->search_type() == SearchTypeRange) {
                RangeSearchParams<T> *p = (RangeSearchParams<T> *) &params;

                if (this->white_range != p->bandwidth || m_index == NULL) {
                    build_index(p->bandwidth);
                }
            }

            // Compile the sample

            typename Point<T>::list *result = new typename Point<T>::list();
            struct kdres *presults;

            // tranform coordinate
            vector<T> xt = this->transform_vector(x);

            // re-package data for KD-Tree
            double pos[this->dimension()];
            for (size_t index = 0; index < this->dimension(); index++) {
                pos[index] = xt[index];
            }

            presults = kd_nearest_range(m_index, pos, this->white_radius);

            // Re-package results
            while (!kd_res_end(presults)) {
                typename Point<T>::ptr ptr = (typename Point<T>::ptr) kd_res_item(presults, &pos[0]);
                result->push_back(ptr);
                kd_res_next(presults);
            }

            // calculate distances
//            if (distances != NULL) {
//                typename Point<T>::list::iterator li;
//                for (li = result->begin(); li != result->end(); li++) {
//                    vector<T> dx = vector_subtract(<#vector<T> &v1#>, <#vector<T> &v2#>);
//                }
//            }

            kd_res_free(presults);

            if (PointIndex<T>::write_index_searches) {
                if (params->search_type() == SearchTypeRange) {
                    RangeSearchParams<T> *p = (RangeSearchParams<T> *) params;
                    this->write_search(x, p->bandwidth, result);
                } else {
                    vector<T> white_ranges(xt.size(), this->white_radius);
                    this->write_search(xt, white_ranges, result);
                }
            }
            return result;
        }
    };
}

#endif
