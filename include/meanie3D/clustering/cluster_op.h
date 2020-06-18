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


#ifndef M3D_CLUSTEROPERATION_H
#define M3D_CLUSTEROPERATION_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/parallel.h>
#include <meanie3D/operations/kernels.h>
#include <meanie3D/clustering/cluster_list.h>
#include <meanie3D/clustering/detection.h>
#include <meanie3D/operations/operation.h>

namespace m3D {

    using namespace boost;

    // Various available post-aggregation methods

    typedef enum
    {
        PostAggregationMethodNone,
        PostAggregationMethodCoalescence,
        PostAggregationMethodDRF
    } PostAggregationMethod;

    /** Encapsulates all steps needed for a complete
     * clustering operation. 
     *
     * TODO: decouple from NetCDF specific code.
     */
    template<typename T>
    class ClusterOperation : public Operation<T>
    {
    private:

        size_t m_cluster_threadcount;
        boost::progress_display *m_progress_bar;
        const detection_params_t<T> m_params;
        const detection_context_t<T> m_context;

    public:

        /** Used to move the progress bar forward if it's switched
         * on during the clustering operation.
         */
        void increment_cluster_progress() const;

        /** Called by threads (if threading is enabled) to report the
         * completion of their separate tasks. When all threads have
         * reported in, the task is complete.
         */
        void report_done();

        ClusterOperation(const detection_params_t<T> &params,
                         const detection_context_t<T> &ctx)
                : Operation<T>(ctx.fs, ctx.index), m_cluster_threadcount(0), m_progress_bar(NULL), m_params(params),
                  m_context(ctx) {};

        virtual ~ClusterOperation() {};

        /** Clustering the whole feature space by running an iterative procedure for each point
         * in feature space.
         *
         * Note, that each point in feature space will be assigned a cluster afterwards. In order
         * to clear this, call reset_clustering() on the feature space.
         * 
         * @param params 
         * 
         * @return ClusterList object with the results.
         */
        ClusterList<T> *cluster();

#pragma mark -
#pragma mark Some public helpers for visualization

#if WRITE_MODES
        static size_t pass_counter();
        static void reset_pass_counter();
        static void increment_pass_counter();
#endif

    };
}

#endif
