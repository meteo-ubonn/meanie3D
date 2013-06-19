#ifndef _M3D_ScalarIndexImpl_H_
#define _M3D_ScalarIndexImpl_H_

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>

namespace m3D {

	using namespace std;
    
    template <typename T>
    ScalarIndex<T>::ScalarIndex(CoordinateSystem<T> *cs)
    : m_coordinate_system(cs)
    , m_data(NULL)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        this->construct_recursive( 0, &m_data, gp );
    }
    
    template <typename T>
    ScalarIndex<T>::ScalarIndex(ScalarIndex<T> *o)
    : m_coordinate_system(o->m_coordinate_system)
    , m_data(NULL)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();

        this->copy_recursive(o,0,gp);
    }
    
    template <typename T>
    ScalarIndex<T>::~ScalarIndex()
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();

        this->destroy_recursive( 0, &m_data, gp );
    }
    
    template <typename T>
    void
    ScalarIndex<T>::copy_recursive(ScalarIndex<T> *otherIndex,
                                   size_t dim_index,
                                   typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
        {
            for ( size_t index = 0; index < dimSize; index++ )
            {
                gridpoint[dim_index] = index;
                
                copy_recursive(otherIndex,dim_index+1,gridpoint);
            }
        }
        else
        {
            typename CoordinateSystem<T>::GridPoint gIter = gridpoint;
            
            for (size_t i=0; i<dimSize; i++)
            {
                gIter[dim_index] = i;
                
                this->set(gIter,otherIndex->get(gIter));
            }
        }
    }

    template <typename T>
    void
    ScalarIndex<T>::construct_recursive(size_t dim_index,
                                       array_t **array,
                                       typename CoordinateSystem<T>::GridPoint &gridpoint,
                                       ScalarIndex<T> *other)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
        {
            // create array
            
            if ( dim_index == 0 )
            {
                vector<void *> *new_array = new vector<void *>(dimSize,NULL);

                *array = new_array;
            
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;
                
                    construct_recursive( dim_index+1, array, gridpoint, other);
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
                    
                    construct_recursive( dim_index+1, &new_array, gridpoint);
                }
            }
        }
        else
        {
            vector<void *> *super_array = *array;

            size_t super_index = gridpoint[dim_index-1];

            vector<T> *new_array = new vector<T>(dimSize,NULL);
            
            super_array->at(super_index) = new_array;
            
            if ( other != NULL )
            {
                typename CoordinateSystem<T>::GridPoint gIter = gridpoint;
                
                for (size_t i=0; i<dimSize; i++)
                {
                    gIter[dim_index] = i;
                    
                    new_array->at(i) = other->get(gIter);
                }
            }
        }
    }
    
    template <typename T>
    void
    ScalarIndex<T>::destroy_recursive(size_t dim_index,
                                      array_t **array,
                                      typename CoordinateSystem<T>::GridPoint &gridpoint)
    {
        NcDim dim = m_coordinate_system->dimensions()[dim_index];
        
        size_t dimSize = dim.getSize();
        
        if (dim_index < (m_coordinate_system->size()-1) )
        {
            if ( dim_index == 0 )
            {
                vector<void *> *the_array = *array;
                
                for ( size_t index = 0; index < dimSize; index++ )
                {
                    gridpoint[dim_index] = index;
                    
                    destroy_recursive( dim_index+1, array, gridpoint);
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

                    destroy_recursive( dim_index+1, &the_array, gridpoint);
                }
                
                delete the_array;
                
                super_array->at(super_index) = NULL;
            }
        }
        else
        {
            vector<void *> *super_array = *array;
            
            size_t super_index = gridpoint[dim_index-1];
            
            vector<T> *the_array = (vector<T> *) super_array->at(super_index);
            
            delete the_array;
            
            super_array->at(super_index) = NULL;
        }
    }
    
#pragma mark -
#pragma mark Accessors
                      
    template <typename T>
    T
    ScalarIndex<T>::get(const typename CoordinateSystem<T>::GridPoint &gp) const
    {
        vector<void *> *array = m_data;
        
        T result = 0;
        
        for (size_t dim_index = 0; dim_index < gp.size(); dim_index++)
        {
            if ( dim_index < gp.size()-1 )
            {
                array = (vector<void *> *) array->at(gp[dim_index]);
            }
            else
            {
                vector<T> *points = (vector<T> *) array;
                
                // Apparently the reverse lookup can hit indexes too
                // high on occasion. This is a problem with the coordinate
                // system class and needs to be fixed eventually
                
                size_t index = gp[dim_index];
                
                if (index < points->size())
                {
                    result = points->at(index);
                }
            }
        }
        
        return result;
    }
    
    template <typename T>
    T
    ScalarIndex<T>::get(const typename CoordinateSystem<T>::Coordinate &c) const
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        m_coordinate_system->reverse_lookup(c,gp);
        
        return get(gp);
    }

    
    template <typename T>
    void
    ScalarIndex<T>::set(const typename CoordinateSystem<T>::GridPoint &gp, T value)
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
                vector<T> *points = (vector<T> *) array;
                
                // Apparently the reverse lookup can hit indexes too
                // high on occasion. This is a problem with the coordinate
                // system class and needs to be fixed eventually
                
                size_t index = gp[dim_index];
                
                if (index < points->size())
                {
                    points->at(index) = value;
                }

            }
        }
    }
    
    template <typename T>
    void
    ScalarIndex<T>::set(const typename CoordinateSystem<T>::Coordinate &c, T value)
    {
        typename CoordinateSystem<T>::GridPoint gp = m_coordinate_system->newGridPoint();
        
        m_coordinate_system->reverse_lookup(c,gp);
        
        set(gp,value);
    }

    
}; //namespace

#endif
