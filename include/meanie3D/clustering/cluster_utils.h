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
            
            /** Iterates over all clusters in the list and replaces the 
             * points in each cluster with the value from the data store. 
             * The method assumes, that the variables from the data store 
             * are in the same order as the variables in the feature 
             * space's value range. This should not be a problem if you 
             * use the data store that the feature space was constructed 
             * from. 
             * 
             * @param list
             * @param dataStore
             */
            static 
            void 
            replace_points_from_datastore(ClusterList<T> &list,
                                          typename DataStore<T>::ptr dataStore);

            /** Iterates over the clusters in the list and checks on each 
             * cluster's points. If a point borders the area marked as
             * off limits in the feature space, the cluster's margin flag
             * is set. 
             * 
             * @param list
             * @param fs
             */
            static
            void obtain_margin_flag(ClusterList<T> &list, 
                                    typename FeatureSpace<T>::ptr fs);
            
    };
}
    
#endif
