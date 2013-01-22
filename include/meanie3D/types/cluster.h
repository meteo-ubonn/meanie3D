#ifndef _M3D_Cluster_H_
#define _M3D_Cluster_H_

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>

#include <vector>
#include <netcdf>

#include <cf-algorithms/cf-algorithms.h>
#include <meanie3D/types/histogram.h>

namespace m3D {
    
    // Forward declarations

	using namespace std;
    using namespace netCDF;
	using cfa::meanshift::Point;
    
    /** Cluster of points in feature space. A cluster is a point in feature space,
     * where multiple trajectories of original feature space points end. This end
     * point is called the cluster's 'mode'.
     */
    template <class T>
    class Cluster
    {
        
    private:
        
        typedef map< size_t, typename Histogram<T>::ptr > histogram_map_t;
        
        histogram_map_t             m_histograms;
        
    public:
        
#pragma mark -
#pragma mark Type definitions / Constants
        
        typedef Cluster<T> *        ptr;
        
        typedef vector<ptr>    		list;
        
        typedef unsigned long long  id_t;
        
        static const id_t           NO_ID;
        
#pragma mark -
#pragma mark Public properties
        
        /** Center of this cluster
         */
        vector<T>   	           	mode;
        
        /** List of feature-space points comprising this cluster
         */
        typename Point<T>::list     points;
        
        /** Unique cluster ID. Used for tracking.
         */
        id_t                        id;
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Default constructor
         */
        Cluster();
        
        /** Constructor
         * @param mode
         */
        Cluster( vector<T> mode );
        
        /** Copy constructor
         * @param cluster
         */
        Cluster( const Cluster<T> &o );
        
        /** Destructor */
        virtual ~Cluster();
        
#pragma mark -
#pragma mark Adding / Removing points
        
        /** Add the feature space point to this cluster
         * @param shared pointer to a feature space point
         */
        void
        add_point( Point<T> *point );
        
        /** Remove the feature space point from this cluster
         * @param feature space point
         */
        void
        remove_point( Point<T> *point );
        
        /** Add a list of points to this cluster.
         * Does not insert duplicates.
         */
        void
        add_points( const vector< Point<T> * > &list );
        
        /** Searches the list of points in this cluster for the given one.
         * @param point
         * @return yes/no
         */
        bool
        has_point( typename Point<T>::ptr point );
        
#pragma mark -
#pragma mark Derived properties
        
        /** Uses the values in the referenced range variable to
         * calculate a center based on the positions of all points
         * weighed by the given variable
         * @param variable index
         * @return centroid
         */
        vector<T>
        weighed_center( const size_t &variable_index );
        
        /** Finds the lower and upper bound of the values in the given variable
         * @param variable index
         * @param lower_bound (return value)
         * @param upper_bound (return value)
         */
        void dynamic_range( const size_t &variable_index, T &lower_bound, T&upper_bound );
        
#pragma mark -
#pragma mark Accessors
        
        size_t size() { return points.size(); };
        
        Point<T> * operator[](size_t index) { return points[index]; };
        
#pragma mark -
#pragma mark Operators
        
        /** Equal operator. Two clusters are considered equal, when they have identical modes.
         * @param cluster
         */
        bool operator == (const Cluster<T> &o);
        
        /** Copy Operator
         * @param cluster
         */
        Cluster<T> operator = (const Cluster<T> &o);

#pragma mark -
#pragma mark Histogram
        
        /** Creates a histogram of the given variable from the points 
         * in this cluster. The histogram is cached. Subsequent calls
         * return the cached histogram, unless clear_histogram_cache()
         * is called first, which will force re-calculation.
         * @param index in point->values to be used
         * @param bottom boundary of variable's values
         * @param top boundary of variable's values
         * @param number of bins in the histogram (default 25)
         * @return handle on the histogram
         */
        const typename Histogram<T>::ptr histogram(size_t variable_index,
                                                   T valid_min,
                                                   T valid_max,
                                                   size_t number_of_bins = 25 );
        
        /** Clears the histogram caches. Subsequent calls to histogram()
         * will return freshly calculated histograms
         */
        void clear_histogram_cache();
        
#pragma mark -
#pragma mark Coverage
        
        /** Counts the number of points in cluster b, that are identical to points
         * in the receiver. Then calculates the percentage of points of in the receiver
         * that are 'covered' in this fashion. The comparison between points is only
         * done using the <b>spatial</b> component.
         *
         * @param another cluster
         * @return coverage in percentage
         */
        float percent_covered_by( const Cluster<T>& b );
    };
    
};
    
#endif
