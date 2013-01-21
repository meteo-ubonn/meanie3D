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
    Cluster<T>::Cluster() : m_histogram(NULL),id(NO_ID)
    {
    }
    
    template <typename T>
    Cluster<T>::Cluster( vector<T> mode ) : m_histogram(NULL), mode(mode), id(NO_ID)
    {
    };
    
    template <typename T>
    Cluster<T>::Cluster( const Cluster<T> &o ) : m_histogram(NULL)
    {
        id = o.id;
        
        mode = o.mode;
        
        points = o.points;
    }
    
    template <typename T>
    Cluster<T>::~Cluster()
    {
        if (m_histogram)
        {
            delete m_histogram;
            
            m_histogram = NULL;
        }
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
            
            M3DPoint<T> *m3p = (M3DPoint<T> *) p;

#if DEBUG_MEANSHIFT_GRAPH
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
    }
    
    template <typename T>
    Cluster<T>
    Cluster<T>::operator = (const Cluster<T> &o)
    {
        return Cluster<T>(o);
    }
    
}; //namespace
    
#endif
