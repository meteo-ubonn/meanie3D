#ifndef _M3D_TRACKING_CLASS_H_
#define _M3D_TRACKING_CLASS_H_

#include <netcdf>

#include <cf-algorithms/id.h>

#include <meanie3D/types/cluster_list.h>
#include <meanie3D/utils.h>

namespace m3D {
    
    using namespace utils;
    using namespace netCDF;
    using cfa::id_t;
    using namespace units::values;
    
    template <typename T>
    class Tracking
    {
    private:
        
        // Member Variables
        
        T   m_dist_weight;  // correlation weight distance
        
        T   m_size_weight;  // correlation weight histogram sum
        
        T   m_corr_weight;  // correlation weight histogram rank correlation
        
        s    m_deltaT;       // What is the time between the given slices (in seconds)
        
        s    m_max_deltaT;   // What is the maximum time between slices for valid tracking (in seconds)

        bool    m_useMeanVelocityConstraint;        // make use of max velocity constraint?
        
        T       m_meanVelocitySecurityPercentage;   // Percentual amount of deviation from mean velocity allowed
        
        meters_per_second   m_maxVelocity;                      // physical maximum speed of objects in m/s
        
        float   m_max_size_deviation;       // how many percent may the objects vary in size (number of points) between scans?

        bool    m_useOverlapConstraint;     // make use of max velocity constraint?
        
        float   m_ms_threshold;             //
        
        float   m_msc_threshold;            // How many percent of coverage is required in splits/merges for track continuation?
        
        
        /** Private default constructor 
         */
        // Tracking() {};
        
    public:
        
        typedef std::vector< Cluster<T> >   track_t;
        
        typedef std::map<id_t, track_t* >   trackmap_t;
        
        /** Constructor
         * @param weight for distance correlation
         * @param weight for size correlation
         * @param weight for histogram rank correlation
         */
        Tracking(T wr=1.0, T wd=1.0, T wt=1.0, s max_delta_t = minute(15))
        : m_dist_weight(wr)
        , m_size_weight(wd)
        , m_corr_weight(wt)
        , m_deltaT(0)   
        , m_max_deltaT(max_delta_t)             // 15 minutes
        , m_useMeanVelocityConstraint(false)    // limit deviation from mean velocity (false)
        , m_meanVelocitySecurityPercentage(0.5) // to 50 %
        , m_maxVelocity(meters_per_second(100.0))                  // limit max velocity to 30 m/s (~108 km/h)
        , m_max_size_deviation(2.5)             // how many percent may the objects vary in size (number of points) between scans (250%)
        , m_useOverlapConstraint(true)
        , m_ms_threshold(0.33)                  // percentage coverage old/new for merge/split (33%)
        , m_msc_threshold(0.75)                 // percentage coverage old/new in merge/split for continuing track (75%)
        {};
        
        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param index of variable to be used for the histogram correlation.
         * @param new clusters, which need id's
         */
        void track(typename ClusterList<T>::ptr previous,
                   typename ClusterList<T>::ptr current,
                   CoordinateSystem<T> *cs,
                   const NcVar &track_variable,
                   Verbosity verbosity = VerbosityNormal);
        
        // Accessors
        
        meters_per_second maxTrackingSpeed() { return m_maxVelocity; }
        void setMaxTrackingSpeed(meters_per_second speed) { m_maxVelocity = speed;}
        
        float mergeSplitThreshold() {return m_ms_threshold;}
        void setMergeSplitThreshold(float value) {m_ms_threshold=value;}
        
        float mergeSplitContinuationThreshold() {return m_msc_threshold;}
        void setMergeSplitContinuationThreshold(float value) {m_msc_threshold=value;}

    private:
        

    };
}

#endif
