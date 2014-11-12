#ifndef M3D_TRACK_H
#define M3D_TRACK_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering.h>

#include <vector>
#include <map>

namespace m3D { 

    /** This simple class constitutes one track. It mainly consists
     * of public properties.
     */        
    template<typename T>
    class Track
    {
        public:
            
#pragma mark -
#pragma mark Type definitions / Constants
            
            typedef Track* ptr;
            typedef std::map< m3D::id_t, ptr> trackmap;

#pragma mark -
#pragma mark Constructor/Destructor

            /** Default constructor.
             */
            Track() {};

#pragma mark -
#pragma mark Public properties

            /** the cluster's identifier
             */
            m3D::id_t id; 

            /** A list of pointers to cluster objects constituting the
             * actual track.
             */
            std::vector< typename Cluster<T>::ptr > clusters; 

            /** list of source files, one per cluster in the list.
             */
            std::vector< std::string > sourcefiles; 

            /** A list of minimum values found for each variable
             * in the cluster's value range.
             */
            std::vector<T> min;

            /** A list of maximum values found for each variable
             * in the cluster's value range.
             */
            std::vector<T> max;
            
#pragma mark -
#pragma mark Public methods
            
            /** Shortcut for clusters.size()
             * 
             * @return number of clusters in this track.
             */
            size_t size() { return clusters.size(); };
            
    };
}

#endif
