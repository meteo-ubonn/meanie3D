#ifndef _M3D_ClusterUtils_H_
#define _M3D_ClusterUtils_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <cf-algorithms/cf-algorithms.h>

#include <meanie3D/types/cluster_list.h>

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
                                          WeightFunction<T> *weight_function,
                                          const Verbosity verbosity = VerbosityNormal);
    };
};
    
#include "cluster_utils_impl.h"
    
#endif
