#ifndef M3D_TRACKING_CLASS_H
#define M3D_TRACKING_CLASS_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/utils.h>
#include <meanie3D/clustering.h>

#include <netcdf>

#include "tracking.h"

namespace m3D {

    using namespace utils;
    
    /** This simple class constitutes one track. It only consists
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
    };

    /** This class contains the tracking code.
     */
    template <typename T>
    class Tracking
    {
        
    private:

        // Member Variables

        T   m_dist_weight;  // correlation weight distance

        T   m_size_weight;  // correlation weight histogram sum

        T   m_corr_weight;  // correlation weight histogram rank correlation

        ::units::values::s    m_deltaT;       // What is the time between the given slices (in seconds)

        ::units::values::s    m_max_deltaT;   // What is the maximum time between slices for valid tracking (in seconds)

        bool    m_useMeanVelocityConstraint;        // make use of max velocity constraint?

        T       m_meanVelocitySecurityPercentage;   // Percentual amount of deviation from mean velocity allowed

        ::units::values::meters_per_second   m_maxVelocity;                      // physical maximum speed of objects in m/s

        float   m_max_size_deviation;       // how many percent may the objects vary in size (number of points) between scans?

        bool    m_useOverlapConstraint;     // make use of max velocity constraint?

        float   m_ms_threshold;             //

        float   m_msc_threshold;            // How many percent of coverage is required in splits/merges for track continuation?


        /** Private default constructor 
         */
        // Tracking() {};

    public:
        
        /** Constructor
         * @param weight for distance correlation
         * @param weight for size correlation
         * @param weight for histogram rank correlation
         */
        Tracking(T wr=1.0, T wd=1.0, T wt=1.0, ::units::values::s max_delta_t = ::units::values::s(930))
        : m_dist_weight(wr)
        , m_size_weight(wd)
        , m_corr_weight(wt)
        , m_deltaT(0)   
        , m_max_deltaT(max_delta_t)             // 15 minutes )(plus 30 seconds slack)
        , m_useMeanVelocityConstraint(false)    // limit deviation from mean velocity (false)
        , m_meanVelocitySecurityPercentage(0.5) // to 50 %
        , m_maxVelocity(::units::values::meters_per_second(100.0))                  // limit max velocity to 30 m/s (~108 km/h)
        , m_max_size_deviation(2.5)             // how many percent may the objects vary in size (number of points) between scans (250%)
        , m_useOverlapConstraint(true)
        , m_ms_threshold(0.33)                  // percentage coverage old/new for merge/split (33%)
        , m_msc_threshold(0.75)                 // percentage coverage old/new in merge/split for continuing track (75%)
        {};

        /** Compares two cluster lists and propagates or assigns new identifiers.
         * @param clusters from the last run
         * @param clusters form the current run
         * @param coordinate system
         * @param name of variable to be used for the histogram correlation. 
         *        If NULL, histogram comparison is ignored completely.
         * @param verbosity
         */
        void track(typename ClusterList<T>::ptr previous,
                   typename ClusterList<T>::ptr current,
                   const CoordinateSystem<T> *cs,
                   const std::string *tracking_variable_name = NULL,
                   Verbosity verbosity = VerbosityNormal);

        // Accessors

        ::units::values::meters_per_second maxTrackingSpeed() { return m_maxVelocity; }
        void setMaxTrackingSpeed(::units::values::meters_per_second speed) { m_maxVelocity = speed;}

        float mergeSplitThreshold() {return m_ms_threshold;}
        void setMergeSplitThreshold(float value) {m_ms_threshold=value;}

        float mergeSplitContinuationThreshold() {return m_msc_threshold;}
        void setMergeSplitContinuationThreshold(float value) {m_msc_threshold=value;}

        ::units::values::s max_deltaT() {return m_max_deltaT;}
        void set_max_deltaT(::units::values::s seconds) {m_max_deltaT = seconds;}

    };
}

#endif
