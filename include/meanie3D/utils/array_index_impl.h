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
    ArrayIndex<T>::ArrayIndex(const FeatureSpace<T> *fs) : m_fs(fs)
    {
        this->construct_array_recursive(0);
    }
    
    template <typename T>
    ArrayIndex<T>::~ArrayIndex()
    {
        this->destroy_array_recursive(0);
    }
    
    template <typename T>
    void
    ArrayIndex<T>::construct_array_recursive(size_t dim_index, void *)
    {
        if ( dim_index < )
    }
    
    template <typename T>
    void
    ArrayIndex<T>::destroy_array_recursive(size_t dim_index)
    {
    }
    
#pragma mark -
#pragma mark Accessors
    
    template <typename T>
    typename Point<T>::ptr
    ArrayIndex<T>::get(const GridPoint &gp)
    {
    }
    
    template <typename T>
    void
    ArrayIndex<T>::set(const GridPoint &gp, typename Point<T>::ptr p)
    {
    }

    
}; //namespace

#endif
