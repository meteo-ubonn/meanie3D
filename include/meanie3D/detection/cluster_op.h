#ifndef _M3D_ClusterOperation_H_
#define _M3D_ClusterOperation_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <boost/progress.hpp>
#include <boost/thread/mutex.hpp>

#include <cf-algorithms/cf-algorithms.h>

#include <meanie3D/types/cluster_list.h>

namespace m3D {

	using namespace cfa::meanshift;
	using namespace boost;

    template <typename T>
    class ClusterOperation : public Operation<T>
    {
        
    private:
        
        boost::progress_display     *m_progress_bar;
        size_t                      m_cluster_threadcount;
        
#ifdef PROVIDE_THREADSAFETY
        mutex                m_cluster_threadcount_mutex;
        mutex                m_clusters_mutex;
        mutex                m_cluster_progress_mutex;
#endif
        
        public :
        
        /** Used to move the progress bar forward if it's switched
         * on during the clustering operation.
         */
        void increment_cluster_progress() const;
        
        /** Called by threads (if threading is enabled) to report the
         * completion of their separate tasks. When all threads have
         * reported in, the task is complete.
         */
        void report_done();
        
        ClusterOperation(FeatureSpace<T> *fs,
        		 	 	 FeatureSpaceIndex<T> *index)
        : Operation<T>( fs, index ), m_progress_bar(NULL) {};
        
        virtual ~ClusterOperation() {}
        
        /** Clustering the whole feature space by running an iterative procedure for each point
         * in feature space.
         *
         * Note, that each point in feature space will be assigned a cluster afterwards. In order
         * to clear this, call reset_clustering() on the feature space.
         *
         * @param params                : range search or knn search parameter
         * @param fuzziness             : absolute, minimum spatial distance between modes.
         * @param kernel                : kernel function for weighing sample points
         * @param weighing function     : externally handed in weighing function, can be NULL
         * @param termcrit_epsilon      : if the distance between iterations gets lower or equal to this value,
         *                                the iteration is stopped.
         * @param termcrit_iterations   : if more than this number of iterations have been done, it stops.
         * @param show_progress_bar     : bool flag. if true, clustering progress is printed out in the console.
         * @return vector of clusters
         */
        ClusterList<T> cluster( const SearchParameters *params,
                               const Kernel<T> *kernel = new GaussianNormalKernel<T>(),
                               const int weight_index = ::cfa::meanshift::MeanshiftOperation<T>::NO_WEIGHT,
                               const double &drf_threshold = 0.95,
                               const bool show_progress_bar = true );
        
#pragma mark -
#pragma mark Some public helpers for visualization

#if WRITE_MODES
        static size_t pass_counter();
        static void reset_pass_counter();
        static void increment_pass_counter();
#endif
        
    };
};
    
#include "cluster_op_impl.h"
    
#endif
