#ifndef _M3D_TRACKING_CLASS_H_
#define _M3D_TRACKING_CLASS_H_

#include <netcdf>

#include <meanie3D/types/cluster_list.h>
#include <meanie3D/utils/verbosity.h>

namespace m3D {
    
    using namespace utils;
    using namespace netCDF;
    
    template <typename T>
    class Tracking
    {
    private:

        // A matrix is a 2D vector construct in this context
        typedef vector< vector<T> > matrix_t;

        /** Creates a new matrix of the given dimensions
         */
        matrix_t create_matrix(size_t width, size_t height);
        
    public:
        
        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param index of variable to be used for the histogram correlation.
         * @param new clusters, which need id's
         */
        void track(const ClusterList<T> &last_clusters,
                   const ClusterList<T> &current_clusters,
                   const vector<NcVar> &feature_variables,
                   const NcVar &track_variable,
                   Verbosity verbosity = VerbosityNormal);
        
    private:
        

    };
}

#endif
