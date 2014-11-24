#ifndef M3D_CLUSTER_H
#define M3D_CLUSTER_H

#include <meanie3D/defines.h>
#include <meanie3D/namespaces.h>
#include <meanie3D/clustering/histogram.h>
#include <meanie3D/units/units.hpp>
#include <meanie3D/featurespace/coordinate_system.h>
#include <meanie3D/weights.h>

#include <vector>
#include <netcdf>

namespace m3D {

    using namespace std;
    using namespace netCDF;

    /** Cluster of points in feature space. A cluster is a point in feature space,
     * where multiple trajectories of original feature space points end. This end
     * point is called the cluster's 'mode'.
     */
    template <class T>
    class Cluster
    {

    private:

        typedef map< size_t, typename Histogram<T>::ptr > histogram_map_t;
        
        typename Point<T>::list     m_points;
        
        bool                        m_has_margin_points;

        histogram_map_t             m_histograms;

        vector<T>                   m_geometrical_center;

        map< size_t, vector<T> >    m_weighed_centers;

        ::units::values::m          m_radius;

        PointIndex<T>               *m_index;

    protected:
        
        size_t                      m_rank;
        size_t                      m_spatial_rank;

        bool                        m_weight_range_calculated;
        T                           m_min_weight;
        T                           m_max_weight;

    public:

#pragma mark -
#pragma mark Type definitions / Constants

        typedef Cluster<T> *    ptr;

        typedef vector<ptr>     list;

#pragma mark -
#pragma mark Public properties

        /** Center of this cluster
         */
        vector<T>                   mode;

        /** Unique cluster ID. Used for tracking.
         */
        m3D::id_t                   id;

        size_t 
        value_rank() { return m_rank - m_spatial_rank; }
        
        size_t
        spatial_rank() { return m_spatial_rank; };

        size_t
        rank() { return m_rank; };

#pragma mark -
#pragma mark Constructor/Destructor

        /** Default Constructor
         */
        Cluster();

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

        /** Returns a reference to the point list.
         * 
         * @return reference to point list
         */
        virtual typename Point<T>::list &get_points();
        
        /** Sets the points given. The existing points
         * are not deleted from memory by default.
         * 
         * @param points
         * @param deletion flag if <code>true</code> the existing points
         * are deleted from memory. If <code>false</code> they are left
         * alone. Make sure you have control over this, as this presents
         * a potentially very large memory leak!
         */
        void set_points(const typename Point<T>::list &points,
                bool delete_existing = false);
        
        /** Returns the number of points in this cluster
         * 
         * @return 
         */
        size_t size() const;

        /** Checks if the cluster has any points 
         * 
         * @return <code>true</code> points list is empty, <code>false</code>
         * otherwise. 
         */
        bool empty() const;
        
        /** Returns the point at the given index.
         * 
         * @param index
         * @return 
         */
        typename Point<T>::ptr operator[](const size_t &index);
        
        /** Returns the point at the given index.
         * 
         * @param index
         * @return 
         */
        typename Point<T>::ptr at(const size_t& index) const;
        
        /** Clears the point list. 
         * 
         * @param deletion_flag if <code>true</code>, the points will be
         * freed from memory. Defaults to <code>false</code> 
         */
        virtual void clear(bool deletion_flag =false);

        /** If true, the cluster has points neighboring points marked
         * as 'off limits'. This means that it's likely that the cluster
         * is only partially visible. This information is important when
         * considering velocity calculations etc.
         * 
         * @return 
         */        
        bool has_margin_points() const;
        
        /** Set the margin flag
         * @param value
         */
        void set_has_margin_points(bool value);

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
         * @return radius in meters. This only returns meters if the
         * spatial dimensions of the point are readily convertible into
         * meters. Otherwise it simply returns the vector length of the
         * average distance vector.
         * @param coordinate system, required for figuring out the correct units.
         */
        ::units::values::m radius(const CoordinateSystem<T> *cs);

#pragma mark -
#pragma mark Value range

        /** Iterates over all points in the cluster and finds
         * the lower and upper bounds of the values inside. 
         * 
         * @param a vector which will contain the minimum values found
         * after the call. One entry per variable in the value range in
         * the same order. 
         * @param a vector which will contain the maximum values found
         * @param a vector which will contain the median values found
         */        
        void 
        variable_ranges(std::vector<T> &min,
                        std::vector<T> &max,
                        std::vector<T> &median);
        
#pragma mark -
#pragma mark DRF merging

        /** Finds the lower and upper bound of the weight function response
         * @param point list
         * @param weight function
         * @param lower_bound (return value)
         * @param upper_bound (return value)
         */
        static void
        dynamic_range(const typename Point<T>::list &list,
                      const WeightFunction<T> *weight_function,
                      T &lower_bound,
                      T &upper_bound );

        /** Finds the lower and upper bound of weight function response in the
         * whole cluster
         * @param cluster
         * @param weight function
         * @param lower_bound (return value)
         * @param upper_bound (return value)
         */
        void
        dynamic_range(const WeightFunction<T> *weight_function,
                      T &lower_bound,
                      T &upper_bound);

#pragma mark -
#pragma mark Coalescence Merging

        /** @return weight function response at the cluster mode
         */
        T modal_weight_response(const WeightFunction<T> *w) const;

        /** @return average weight function response on all points of the
         * cluster
         */
        T average_weight_response(const WeightFunction<T> *w) const;
    };
}
    
#endif
