#ifndef _M3D_TRACKING_CLASS_H_
#define _M3D_TRACKING_CLASS_H_

#include <meanie3D/types/cluster_list.h>

namespace m3D {
    
    template <typename T>
    class Tracking
    {
        
    public:
        
        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param new clusters, which need id's
         */
        void track(const ClusterList<T> &last_clusters,
                   const ClusterList<T> &current_clusters );
    };
}

#endif
