#ifndef _M3D_Cluster_Impl_H_
#define _M3D_Cluster_Impl_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <meanie3D/types/point.h>

namespace m3D {

	using namespace std;
	using cfa::meanshift::Point;

#pragma mark -
#pragma mark Constructor/Destructor
    
    template <typename T>
    Cluster<T>::Cluster()
    : m_radius(-1)
    , m_index(NULL)
    , m_dimension(0)
    , m_spatial_dimension(0)
    , m_weight_range_calculated(false)
    , m_min_weight(0)
    , m_max_weight(0)
    , id(NO_ID)
    {}
    
    template <typename T>
    Cluster<T>::Cluster(const vector<T> &mode, size_t spatial_dimension)
    : m_radius(-1)
    , m_index(NULL)
    , m_dimension(mode.size())
    , m_spatial_dimension(spatial_dimension)
    , m_weight_range_calculated(false)
    , m_min_weight(0)
    , m_max_weight(0)
    , mode(mode)
    , id(NO_ID)
    {
        assert( m_dimension > m_spatial_dimension ); // TODO: change on refac #146
    }
    
    template <typename T>
    Cluster<T>::Cluster( const Cluster<T> &o )
    : m_radius(-1)
    , m_index(NULL)
    , m_dimension(o.m_dimension)
    , m_spatial_dimension( o.m_spatial_dimension )
    , m_weight_range_calculated(o.m_weight_range_calculated)
    , m_min_weight(o.m_min_weight)
    , m_max_weight(o.m_max_weight)
    , mode(o.mode)
    , points(o.points)
    , id(o.id)
    {
    }
    
    template <typename T>
    Cluster<T>::~Cluster()
    {
        this->clear_histogram_cache();

        this->clear_index();
    }
    
#pragma mark -
#pragma mark Adding / Removing points
    
    template <typename T>
    void
    Cluster<T>::add_point( Point<T> *point )
    {
        points.push_back( point );
        
#if PROVIDE_THREADSAFETY
        point.mutex.lock();
#endif
        M3DPoint<T> *m3p = (M3DPoint<T> *) point;

        m3p->cluster = this;
        
#if PROVIDE_THREADSAFETY
        point.mutex.unlock();
#endif
    }
    
    template <typename T>
    void
    Cluster<T>::add_points( const vector< Point<T> * > &list, bool addOriginalPointsOnly )
    {
        // cout << "Cluster::add_points() " << this->mode << " -> #points = " << this->points.size() << endl;
        
        typename Point<T>::list::const_iterator li;
        
        for ( li = list.begin(); li != list.end(); li++ )
        {
            typename Point<T>::ptr p = *li;
            
            M3DPoint<T> *m3p = (M3DPoint<T> *) p;

#if DEBUG_MEANSHIFT_GRAPH
            if ( m3p->cluster )
            {
                cerr << "Point " << m3p->values << " is already assigned to cluster at " << m3p->cluster->mode << endl;
                continue;
            }
#endif
            if ( !m3p->isOriginalPoint && addOriginalPointsOnly )
            {
                continue;
            }
            
            this->add_point( p );
        }
    }
    
    template <typename T>
    void
    Cluster<T>::remove_point( Point<T> *point )
    {
        typename vector< Point<T> * >::iterator it = points.find( *point );
        
        if ( it != points.end() )
        {
            points.erase( it );
            
#if PROVIDE_THREADSAFETY
            point->mutex.lock();
#endif
            point->cluster = NULL;
            
#if PROVIDE_THREADSAFETY
            point->mutex.unlock();
#endif
            
        }
    }
    
    template <typename T>
    bool
    Cluster<T>::has_point( typename Point<T>::ptr point )
    {
        typename Point<T>::list::iterator fi = find( points.begin(), points.end(), point );
        
        return fi != points.end();
    }
    
#pragma mark -
#pragma mark Derived properties
    
    template <typename T>
    vector<T>
    Cluster<T>::weighed_center( const size_t &variable_index )
    {
        throw "not implemented";
    }
    
#pragma mark -
#pragma mark Operators
    
    template <typename T>
    bool
    Cluster<T>::operator == (const Cluster<T> &o)
    {
        return o.mode() == mode();
    }
    
    template <typename T>
    Cluster<T>
    Cluster<T>::operator = (const Cluster<T> &o)
    {
        return Cluster<T>(o);
    }
    
#pragma mark -
#pragma mark Histograms
    
    template <typename T>
    const typename Histogram<T>::ptr
    Cluster<T>::histogram(size_t variable_index,
                          T valid_min,
                          T valid_max,
                          size_t number_of_bins)
    {
        Histogram<T> *h = NULL;
        
        try
        {
            h = this->m_histograms.at(variable_index);
        }
        catch (const std::exception& e)
        {
            h = Histogram<T>::create( this->points, variable_index, valid_min, valid_max, number_of_bins );
            
            this->m_histograms.insert( std::pair< size_t, typename Histogram<T>::ptr >( variable_index, h ) );
        }
        
        return h;
    }
    
    template <typename T>
    void
    Cluster<T>::clear_histogram_cache()
    {
        typename histogram_map_t::iterator it;
        
        for ( it = this->m_histograms.begin(); it != this->m_histograms.end(); ++it )
        {
            delete it->second;
        }
        
        this->m_histograms.clear();
    }
    
#pragma mark -
#pragma mark Index
    
    template <typename T>
    PointIndex<T> *
    Cluster<T>::index()
    {
        if ( this->m_index == NULL )
        {
            // index is constructed using the spatial components only
            
            vector<size_t> indexes(m_spatial_dimension);
            
            for (size_t i=0; i<m_spatial_dimension; i++) indexes[i]=i;

            this->m_index = PointIndex<T>::create( &this->points, indexes );
        }
        
        return m_index;
    }
    
    template <typename T>
    void
    Cluster<T>::clear_index()
    {
        if ( this->m_index != NULL )
        {
            delete this->m_index;
            
            this->m_index = NULL;
        }
    }

#pragma mark -
#pragma mark Coverage
    
    template <typename T>
    vector<T>
    Cluster<T>::geometrical_center(size_t spatial_dimensions)
    {
        if ( m_geometrical_center.empty() )
        {
            vector<T> center(spatial_dimensions,0.0);
        
            typename Point<T>::list::iterator pi;
            
            for ( pi = this->points.begin(); pi != this->points.end(); ++pi )
            {
                vector<T> tmp = center;
                
                center += (*pi)->coordinate;
                
                if ( isnan( center[0] ) )
                {
                    cerr << "ARGHH!" << endl;
                }
            }
            
            m_geometrical_center = center/((T)this->points.size());
        }
        
        return m_geometrical_center;
    }
    
    template <typename T>
    vector<T>
    Cluster<T>::weighed_center(size_t spatial_dimensions, size_t variable_index)
    {
        vector<T> wc;
        
        try
        {
            wc = this->m_weighed_centers.at(variable_index);
        }
        catch (const std::exception& e)
        {
            // first, figure out the min/max values
            
            typename Point<T>::list::const_iterator pi;
            
            T min = std::numeric_limits<T>::max();
            T max = std::numeric_limits<T>::min();
            for ( pi = this->points.begin(); pi != this->points.end(); pi++ )
            {
                T val = (*pi)->values[variable_index];
                if (val < min) {
                    min = val;
                }
                if (val>max) {
                    max = val;
                }
            }
            
            T range = max - min;
            
            // weight is assigned based on percentage of the range
            // within [min..max]
            
            vector<T> center(spatial_dimensions,0.0);
            
            T overall_mass_percent = 0.0;
            
            for ( pi = this->points.begin(); pi != this->points.end(); pi++ )
            {
                T mass = (*pi)->values[variable_index];
                
                T massPercent = (max-mass)/range;
                
                center += (*pi)->coordinate * massPercent;
                
                overall_mass_percent += massPercent;
            }
            
            center /= overall_mass_percent;
            
            m_weighed_centers.insert( std::pair< size_t, vector<T> >( variable_index, center ) );
            
            wc = m_weighed_centers[variable_index];
        }
        
        return wc;
    }
    
    template <typename T>
    void
    Cluster<T>::clear_center_caches()
    {
        m_geometrical_center.clear();
        
        m_weighed_centers.clear();
    }
    
    template <typename T>
    ::units::values::m
    Cluster<T>::radius(const CoordinateSystem<T> *cs)
    {
        if ( m_radius == ::units::values::m(-1) )
        {
            vector<T> mode_spatial( &mode[0], &mode[m_spatial_dimension] );
            
            vector<T> mode_in_meters = cs->to_meters(mode_spatial);
            
            typename Point<T>::list::iterator pi;
            
            ::units::values::m sum;

            for ( pi = this->points.begin(); pi != this->points.end(); pi++ )
            {
                typename Point<T>::ptr p = *pi;
                
                vector<T> coord_in_meters = cs->to_meters(p->coordinate);
                
                ::units::values::m dist = ::units::values::m(vector_norm( coord_in_meters - mode_in_meters ));
                
                sum += dist;
            }
            
            m_radius = sum / ((T)this->points.size());
        }
        
        return m_radius;
    }
    
    template <typename T>
    T
    Cluster<T>::modal_weight_response(const WeightFunction<T> *w) const
    {
        T result = 0;
        
        // pick the first point of this cluster to figure out the
        // spatial and value dimensions
        
        if (!points.empty() && w!=NULL)
        {
            Point<T> *p = points[0];
            
            size_t spatial_dim = p->coordinate.size();
            
            vector<T> mode_spatial_coord( &mode[0], &mode[ spatial_dim ] );
            
            result = w->operator()(mode_spatial_coord);
        }
        
        return result;
    }
    
    template <typename T>
    T
    Cluster<T>::average_weight_response(const WeightFunction<T> *w) const
    {
        T result = 0;
        
        // pick the first point of this cluster to figure out the
        // spatial and value dimensions
        
        for (size_t i = 0; i < points.size(); i++)
        {
            Point<T> *p = points[i];
            
            result += w->operator()(p->coordinate);
        }
        
        if (!points.empty())
        {
            result /= ((T)points.size());
        }
    
        return result;
    }
    
#pragma mark -
#pragma mark Dynamic Range Calculation
    
    template <class T>
    void
    Cluster<T>::dynamic_range(const typename Point<T>::list &points,
                              const WeightFunction<T> *weight_function,
                              T &lower_bound,
                              T&upper_bound )
    {
        lower_bound = std::numeric_limits<T>::max();
        
        upper_bound = std::numeric_limits<T>::min();
        
        typename Point<T>::list::const_iterator pi;
        
        for ( pi = points.begin(); pi != points.end(); pi++ )
        {
            typename Point<T>::ptr p = *pi;
            
            T value = weight_function->operator()(p);
            
            if ( value < lower_bound )
            {
                lower_bound = value;
            }
            
            if ( value > upper_bound )
            {
                upper_bound = value;
            }
        }
    }
    
    template <typename T>
    void
    Cluster<T>::dynamic_range(const WeightFunction<T> *weight_function,
                              T &lower_bound,
                              T &upper_bound)
    {
        if (!m_weight_range_calculated)
        {
            Cluster<T>::dynamic_range( this->points, weight_function, m_min_weight, m_max_weight );
            m_weight_range_calculated=true;
        }
        
        lower_bound = m_min_weight;
        upper_bound = m_max_weight;
    }



} //namespace

#endif
