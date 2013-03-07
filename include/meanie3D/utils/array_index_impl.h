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
    ArrayIndex<T>::construct_array_recursive(size_t dim_index, void **data, typename CoordinateSystem<T>::GridPoint &gridpoint )
    {
        NcDim dim = m_fs->coordinate_system()->dimensions()[dim_index];
        
        if (dim_index==0)
        {
            *data = (void *) new vector<void *>(dim.getSize(),NULL);
            
            gridpoint[

            construct_array_recursive(dim_index+1, data, gridpoint);

        }
        else
        {
            for ( int index = 0; index < dim.getSize(); index++ )
            {
                gridpoint[dimensionIndex] = index;

                vector<void *> *array = (vector<void *>) (*data);

                if (dim_index < (m_fs->coordinate_system()->size()-1) )
                {
                    else
                    {
                        array[index] = new vector<void *>(dim.getSize(),NULL);
                    }
                    
                    construct_array_recursive(dim_index+1, data, gridpoint);
                }
                else
                {
                    array[index] = new vector<typename Point<T>::ptr>(dim.getSize(),NULL);
                }
            }
        }
    }
    
    template <typename T>
    void
    ArrayIndex<T>::destroy_array_recursive(size_t dim_index,typename CoordinateSystem<T>::GridPoint &gridpoint)
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
