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

#ifndef M3D_CLUSTERINDEX_IMPL_H
#define M3D_CLUSTERINDEX_IMPL_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/featurespace/point.h>

#include "cluster_index.h"

namespace m3D {
    namespace utils {

        template <typename T>
        ClusterIndex<T>::ClusterIndex(typename Cluster<T>::list &list,
                const vector<size_t> &dimensions)
        {
            // Figure out the rank of the multiarray needed
            // from the first cluster's mode

            this->m_index = new MultiArrayBlitz<m3D::id_t>(dimensions, m3D::NO_ID);

            if (!list.empty()) {
                for (size_t ci = 0; ci < list.size(); ci++) {
                    typename Cluster<T>::ptr c = list.at(ci);
                    for (size_t pi = 0; pi < c->size(); pi++) {
                        typename Point<T>::ptr p = c->at(pi);
                        this->m_index->set(p->gridpoint, c->id);
                    }
                }
            }
        }

        template <typename T>
        ClusterIndex<T>::~ClusterIndex()
        {
            if (m_index != NULL) {
                delete m_index;
                m_index = NULL;
            }
        }

        template <typename T>
        typename ClusterIndex<T>::index_t *
        ClusterIndex<T>::data()
        {
            return this->m_index;
        }

        template <typename T>
        size_t
        ClusterIndex<T>::count_common_points(const ClusterIndex<T> *a,
                const ClusterIndex<T> *b,
                ::m3D::id_t id)
        {
            assert(a->data()->rank() == b->data()->rank());
            size_t common_points;
            size_t dim = a->data().rank();
            vector<int> gridpoint(dim);
            if (dim == 2) {
                for (size_t i1 = 0; i1 < a->data()->get_dimensions(0); i1++) {
                    gridpoint[0] = i1;
                    for (size_t i2 = 0; i2 < a->data()->get_dimensions(1); i2++) {
                        gridpoint[1] = i2;
                        m3D::id_t id_a = a->data()->get(gridpoint);
                        if (id_a == id) {
                            m3D::id_t id_b = b->data()->get(gridpoint);
                            if (id_b == id_a) common_points++;
                        }
                    }
                }
            } else {
                throw std::runtime_error("Not implemented");
            }
            return common_points;
        }

        template <typename T>
        double
        ClusterIndex<T>::occupation_ratio(const Cluster<T> *cluster_a,
                const Cluster<T> *cluster_b) const
        {
            int common_points = 0;
            for (int i = 0; i < cluster_a->size(); i++) {
                typename Point<T>::ptr p = cluster_a->at(i);
                if (m_index->get(p->gridpoint) == cluster_b->id) {
                    common_points++;
                }
            }
            double ratio = ((double) common_points) / ((double) cluster_a->size());
            return ratio;
        }
    }
}

#endif
