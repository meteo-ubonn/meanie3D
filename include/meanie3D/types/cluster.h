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
    using cfa::meanshift::PointIndex;

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
        
        vector<T>                   m_geometrical_center;
        
        map< size_t, vector<T> >    m_weighed_centers;
        
        T                           m_radius;
        
        PointIndex<T>               *m_index;
        
        size_t                      m_dimension;                // dimension of the points
        
        size_t                      m_spatial_dimension;        // number of spatial dimensions (always first)
        
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
        
        // default constructor
        Cluster();
        
        size_t
        spatial_dimension() { return m_spatial_dimension; };
        
        size_t
        dimension() { return m_dimension; };
        
#pragma mark -
#pragma mark Constructor/Destructor
        
        /** Constructor
         * @param the cluster mode in feature-space (not in spatial coordinates!!)
         * @param number of spatial dimensions
         */
        Cluster(const vector<T> &mode, size_t spatial_dimension);
        
        /** Copy constructor
         * @param cluster
         */
        Cluster(const Cluster<T> &o);
        
        /** Destructor 
         */
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
        add_points( const vector< Point<T> * > &list, bool addOriginalPointsOnly=true );
        
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
#pragma mark Index

        /** Returns the cluster's index (an index over the cluster's points).
         * The result is cached. The cache can be cleared by calling clear_index().
         * @return index on the cluster's points.
         */
        PointIndex<T> *index();

        /** Deletes the cluster's index. Subsequent call to index() will force
         * a re-calculation of the index.
         */
        void clear_index();

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
        float percent_covered_by( const Cluster<T>::ptr b ) const;
        
#pragma mark -
#pragma mark Derived Properties
        
        /** Calculate the cluster center. 
         * The result is cached. The cached value can be reset by
         * calling clear_center_caches();
         * @param number of spatial dimensions
         * @return the spatial coordinate of the geometrical center 
         * (using only the spatial components for calculation)
         */
        vector<T> geometrical_center(size_t spatial_dimensions);
        
        /** Calculate the cluster center weighed by a variable. 
         * The result is cached. The cached value can be reset by 
         * calling clear_center_caches();
         * @param number of spatial dimensions
         * @param index of the variable used for weighing
         * @return the spatial coordinate of the geometrical center
         * (using only the spatial components for calculation)
         */
        vector<T> weighed_center(size_t spatial_dimensions, size_t variable_index);
        
        /** Erases the calculation results of all center calculations,
         * forcing a fresh calculation the next time they are called.
         */
        void clear_center_caches();
        
        /** Returns the maximum distance between the cluster mode and
         * the point farthest away from it.
         */
        T radius();
    };
    
};
    
#endif
