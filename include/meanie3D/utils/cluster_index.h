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

#ifndef M3D_CLUSTERINDEX_H
#define M3D_CLUSTERINDEX_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/array.h>
#include <meanie3D/clustering.h>

namespace m3D { namespace utils {

    template <typename T>
    struct ClusterIndex
    {

    public:

        typedef MultiArray<m3D::id_t> index_t;

    private:

        index_t *m_index;

    public:

        ClusterIndex(typename Cluster<T>::list &list,
                     const vector<size_t> &dimensions);

        ~ClusterIndex();

        // Accessors

        index_t *data();

        // Static functions

        /** Iterates over the other array, finds all point
         */
        static size_t count_common_points(const ClusterIndex<T> *a,
                                          const ClusterIndex<T> *b,
                                          ::m3D::id_t id);

        /** What is the ratio of gridpoints in cluster A that
         * are also occupied by points from cluster B? 
         *
         * Important: cluster B must have been indexed by this 
         * index prior to calling.
         *
         * @param cluster A
         * @param cluster B
         */
        double occupation_ratio(const Cluster<T> *cluster_a,
                                const Cluster<T> *cluster_b) const;

    };
}}

#endif
