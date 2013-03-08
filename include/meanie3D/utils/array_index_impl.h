#ifndef _M3D_ArrayIndex_Impl_H_
#define _M3D_ArrayIndex_Impl_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

#include <meanie3D/types/point.h>

namespace m3D {

	using namespace std;
	using cfa::meanshift::Point;
    

    template <typename T>
    ArrayIndex<T>::ArrayIndex(const FeatureSpace<T> *fs) : m_fs(fs), m_data(NULL)
    {
        typename CoordinateSystem<T>::GridPoint gp = fs->coordinate_system->newGridPoint();
        
        this->construct_array_recursive( 0, &m_data, gp );
        
        this->build_index();
    }
    
    template <typename T>
    ArrayIndex<T>::~ArrayIndex()
    {
        typename CoordinateSystem<T>::GridPoint gp = m_fs->coordinate_system->newGridPoint();

        this->destroy_array_recursive( 0, &m_data, gp );
    }
    
    template <typename T>
    void
    ArrayIndex<T>::construct_array_recursive(size_t dim_index, array_t **array, typename CoordinateSystem<T>::GridPoint &gridpoint )
    {
        NcDim dim = m_fs->coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_fs->coordinate_system->size()-1) )
        {
            // create array
            
            if ( dim_index == 0 )
            {
                vector<void *> *new_array = new vector<void *>(dimSize,NULL);

                *array = new_array;
            
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;
                
                    construct_array_recursive( dim_index+1, array, gridpoint);
                }
            }
            else
            {
                vector<void *> *super_array = *array;
                
                size_t super_index = gridpoint[dim_index-1];
                
                vector<void *> *new_array = new vector<void *>(dimSize,NULL);
                
                super_array->at(super_index) = new_array;

                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;
                    
                    construct_array_recursive( dim_index+1, &new_array, gridpoint);
                }
            }
        }
        else
        {
            vector<void *> *super_array = *array;

            size_t super_index = gridpoint[dim_index-1];

            vector<typename Point<T>::ptr> *new_array = new vector<typename Point<T>::ptr>(dimSize,NULL);
            
            super_array->at(super_index) = new_array;
        }
    }
    
    template <typename T>
    void
    ArrayIndex<T>::destroy_array_recursive(size_t dim_index, array_t **array, typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = m_fs->coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_fs->coordinate_system->size()-1) )
        {
            if ( dim_index == 0 )
            {
                vector<void *> *the_array = *array;
                
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;
                    
                    destroy_array_recursive( dim_index+1, array, gridpoint);
                }
                
                delete the_array;
                
                *array = NULL;
            }
            else
            {
                size_t super_index = gridpoint[dim_index-1];
                
                vector<void *> *super_array = *array;
                
                vector<void *> *the_array = (vector<void *> *) super_array->at(super_index);
                
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;

                    destroy_array_recursive( dim_index+1, &the_array, gridpoint);
                }
                
                delete the_array;
                
                super_array->at(super_index) = NULL;
            }
        }
        else
        {
            vector<void *> *super_array = *array;
            
            size_t super_index = gridpoint[dim_index-1];
            
            vector<typename Point<T>::ptr> *the_array = (vector<typename Point<T>::ptr> *) super_array->at(super_index);
            
            delete the_array;
            
            super_array->at(super_index) = NULL;
        }
    }
    
    template <typename T>
    void
    ArrayIndex<T>::build_index()
    {
        CoordinateSystem<T> *cs = this->m_fs->coordinate_system;
        
        for (size_t i=0; i<m_fs->points.size(); i++)
        {
            typename Point<T>::ptr p = m_fs->points[i];
            
            this->set( p->gridpoint, p );
        }
    }

    
#pragma mark -
#pragma mark Accessors
                      
    template <typename T>
    typename Point<T>::ptr
    ArrayIndex<T>::get(const typename CoordinateSystem<T>::GridPoint &gp)
    {
        vector<void *> *array = m_data;
        
        typename Point<T>::ptr result = NULL;
        
        for (size_t dim_index = 0; dim_index < gp.size(); dim_index++)
        {
            if ( dim_index < gp.size()-1 )
            {
                array = array->at(gp[dim_index]);
            }
            else
            {
                vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;
                
                result = points->at(gp[dim_index]);
            }
        }
        
        return result;
    }
    
    template <typename T>
    void
    ArrayIndex<T>::set(const typename CoordinateSystem<T>::GridPoint &gp, typename Point<T>::ptr p)
    {
        vector<void *> *array = m_data;
        
        for (size_t dim_index = 0; dim_index < gp.size(); dim_index++)
        {
            if ( dim_index < gp.size()-1 )
            {
                array = (vector<void *> *) array->at(gp[dim_index]);
            }
            else
            {
                vector<typename Point<T>::ptr> *points = (vector<typename Point<T>::ptr> *) array;
                
                points->at(gp[dim_index]) = p;
            }
        }
    }

    
}; //namespace

#endif
