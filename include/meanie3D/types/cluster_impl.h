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

    template <typename T>
    const typename Cluster<T>::id_t Cluster<T>::NO_ID = numeric_limits<id_t>::max();
    
#pragma mark -
#pragma mark Constructor/Destructor
    
    template <typename T>
    Cluster<T>::Cluster() : id(NO_ID)
    {
    }
    
    template <typename T>
    Cluster<T>::Cluster( vector<T> mode ) : mode(mode), id(NO_ID)
    {
    };
    
    template <typename T>
    Cluster<T>::Cluster( const Cluster<T> &o )
    {
        id = o.id;
        
        mode = o.mode;
        
        points = o.points;
    }
    
    template <typename T>
    Cluster<T>::~Cluster()
    {
        this->clear_histogram_cache();
    };
    
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
    Cluster<T>::add_points( const vector< Point<T> * > &list )
    {
        // cout << "Cluster::add_points() " << this->mode << " -> #points = " << this->points.size() << endl;
        
        typename Point<T>::list::const_iterator li;
        
        for ( li = list.begin(); li != list.end(); li++ )
        {
            typename Point<T>::ptr p = *li;
            
#if DEBUG_MEANSHIFT_GRAPH
            M3DPoint<T> *m3p = (M3DPoint<T> *) p;

            if ( m3p->cluster )
            {
                cerr << "Point " << m3p->values << " is already assigned to cluster at " << m3p->cluster->mode << endl;
                continue;
            }
#endif
            typename Point<T>::list::iterator fi = find( points.begin(), points.end(), p );
            
            if ( fi == points.end() )
            {
                this->add_point( p );
            }
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
    };
    
    template <typename T>
    bool
    Cluster<T>::has_point( typename Point<T>::ptr point )
    {
        typename Point<T>::list::iterator fi = find( points.begin(), points.end(), point );
        
        return fi != points.end();
    };
    
#pragma mark -
#pragma mark Derived properties
    
    template <typename T>
    vector<T>
    Cluster<T>::weighed_center( const size_t &variable_index )
    {
        throw "not implemented";
    }
    
    template <typename T>
    void
    Cluster<T>::dynamic_range( const size_t &variable_index, T &lower_bound, T&upper_bound )
    {
        return Point<T>::dynamic_range( points, variable_index, lower_bound, upper_bound );
    }
    
#pragma mark -
#pragma mark Operators
    
    template <typename T>
    bool
    Cluster<T>::operator == (const Cluster<T> &o)
    {
        return o.mode() == mode();
    };
    
    template <typename T>
    Cluster<T>
    Cluster<T>::operator = (const Cluster<T> &o)
    {
        return Cluster<T>(o);
    };
    
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
#pragma mark Coverage
    
    template <typename T>
    float
    Cluster<T>::percent_covered_by( const Cluster<T>& b )
    {
        // TODO: this can be done a lot faster using an index
        
        size_t num_common_points = 0;
        
        typename Point<T>::list::iterator a_points;
        
        for ( a_points = this->points.begin(); a_points != this->points.end(); a_points++ )
        {
            typename Point<T>::ptr pa = *a_points;

            typename Point<T>::list::iterator b_points;

            for ( b_points = b.points.begin(); b_points != b.points.end(); b_points++ )
            {
                typename Point<T>::ptr pb = *b_points;
                
                if ( pa->coordinate == pb->coordinate )
                {
                    num_common_points++;
                }
            }
        }
        
        return ((float)num_common_points) / ((float)this->points.size());
    }
    
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
            vector<T> center(spatial_dimensions,0.0);
            
            typename Point<T>::list::iterator pi;
            
            T overall_mass = 0.0;
            
            for ( pi = this->points.begin(); pi != this->points.end(); pi++ )
            {
                T mass = (*pi)->values[variable_index];
                
                center += (*pi)->coordinate * mass;
                
                overall_mass += mass;
            }
            
            center /= overall_mass;
            
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

    
    
}; //namespace

#endif
