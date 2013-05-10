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
        
        // Matrix

        // A matrix is a 2D vector construct in this context
        typedef vector< vector<T> > matrix_t;
        
        /** Creates a new matrix of the given dimensions
         */
        matrix_t
        create_matrix(size_t width, size_t height);
        
        // Flag Matrix

        typedef vector< vector<int> > flag_matrix_t;
        
        /** Creates a new matrix of the given dimensions
         */
        flag_matrix_t
        create_flag_matrix(size_t width, size_t height, int defaultValue=0 );
        
        
        // Member Variables
        
        
        T   m_dist_weight;  // correlation weight distance
        
        T   m_size_weight;  // correlation weight histogram sum
        
        T   m_corr_weight;  // correlation weight histogram rank correlation
        
        T   m_deltaT;       // This should in future be calculated properly from timestamps

        bool    m_useMeanVelocityConstraint;        // make use of max velocity constraint?
        
        T       m_meanVelocitySecurityPercentage;   // Percentual amount of deviation from mean velocity allowed
        
        T       m_maxVelocity;                      // physical maximum speed of objects in m/s
        
        bool    m_useOverlapConstraint;             // make use of max velocity constraint?
        
        float   m_merge_threshold;
        
        float   m_max_size_deviation;               // how many percent may the objects vary in size (number of points) between scans?
        
        
        /** Private default constructor 
         */
        // Tracking() {};
        
    public:
        
        typedef std::vector< Cluster<T> >                   track_t;
        
        typedef std::map< cfa::meanshift::id_t, track_t* >  trackmap_t;
        
        /** Constructor
         * @param weight for distance correlation
         * @param weight for size correlation
         * @param weight for histogram rank correlation
         */
        Tracking(T wr=1.0, T wd=1.0, T wt=1.0 )
        : m_dist_weight(wr)
        , m_size_weight(wd)
        , m_corr_weight(wt)
        , m_deltaT(300)                         // 5 minutes
        , m_useMeanVelocityConstraint(false)    // limit deviation from mean velocity (false)
        , m_meanVelocitySecurityPercentage(0.5) // to 50 %
        , m_maxVelocity(100.0)                  // limit max velocity to 30 m/s (~108 km/h)
        , m_useOverlapConstraint(true)
        , m_merge_threshold(0.33)               // percentage coverage old/new for merge/split (33%)
        , m_max_size_deviation(2.5)             // how many percent may the objects vary in size (number of points) between scans (250%)
        {};
        
        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param index of variable to be used for the histogram correlation.
         * @param new clusters, which need id's
         */
        void track(typename ClusterList<T>::ptr previous,
                   typename ClusterList<T>::ptr current,
                   const NcVar &track_variable,
                   Verbosity verbosity = VerbosityNormal);
        
    private:
        

    };
}

#endif
