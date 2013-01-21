#ifndef _M3D_TRACKING_CLASS_Impl_H_
#define _M3D_TRACKING_CLASS_Impl_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>

#include <meanie3D/types/cluster.h>
#include <meanie3D/tracking/tracking.h>


namespace m3D {

    template <typename T>
    void
    Tracking<T>::track(const ClusterList<T> &last_clusters,
                       const ClusterList<T> &current_clusters )
    {
        size_t current_id = Cluster<T>::NO_ID;

        // Find out the highest id from the previous list
        
        typename Cluster<T>::list::iterator ci;
        
        typename Cluster<T>::list previous = last_clusters.clusters;
        
        for ( ci = previous.begin(); ci != previous.end(); ci++ )
        {
            typename Cluster<T>::ptr c = *ci;
            
            if ( c->id > current_id )
            {
                current_id = c->id;
            }
        }
        
        current_id++;
        
#if DEBUG_TRACKING
        cout << "Tracking: next available id is " << current_id << endl;
#endif

    }
    
};

#endif
