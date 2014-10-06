#ifndef M3D_CLUSTERUTILS_H
#define M3D_CLUSTERUTILS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <meanie3D/clustering/cluster.h>
#include <meanie3D/clustering/cluster_list.h>

namespace m3D {

    template <typename T>
    class ClusterUtils
    {
        private:

            float   m_merge_threshold;

        public:

            ClusterUtils(float merge_threshold = 0.33);

            /** Calculates which of the new clusters are the result of 
             * splitting of two previous clusters and merges them again.
             *
             * @param clusters from the last run
             * @param clusters from the current run
             * @param weight function
             */
            void
            filter_with_previous_clusters(typename ClusterList<T>::ptr previous,
                                          typename ClusterList<T>::ptr current,
                                          CoordinateSystem<T> *coord_system,
                                          WeightFunction<T> *weight_function,
                                          const Verbosity verbosity = VerbosityNormal);
    };
}
    
#endif
